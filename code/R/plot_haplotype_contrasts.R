#!/usr/bin/env Rscript
# Native R participant-carrier contrasts. Observed differences and inference are separate.
# Missing/unavailable values remain NA; no test is run by this renderer.
entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
if (!exists("read_tsv", mode = "function")) source(file.path(dirname(normalizePath(entry)), "compare_haplotypes.R"))

CONTRAST_FIELDS <- strsplit("locus_id chrom lead_pos haplotype_id contrast case_group control_group panel_status comparison_status inference_status reason n_resolved_case n_resolved_control carriers_case carriers_control carrier_frequency_case carrier_frequency_control descriptive_carrier_frequency_difference carrier_frequency_difference CI95_low CI95_high p_fisher_two_sided q_BH_all_pairwise_tests haplotype_frequency_case_descriptive haplotype_frequency_control_descriptive", " ")[[1]]
CONTRAST_SOURCES <- c("run_summary.json", "lead_loci.tsv", "haplotype_comparisons.tsv", "haplotype_group_frequencies.tsv", "COMPLETED.txt")
CONTRAST_INFERENCE_FIELDS <- c("conditional_odds_ratio_a_vs_b", "CI95_low", "CI95_high", "p_fisher_two_sided", "q_BH_all_haplotype_tests", "carrier_frequency_difference_a_minus_b", "carrier_difference_CI95_low", "carrier_difference_CI95_high", "q_BH_all_pairwise_tests")

contrast_numeric <- function(row, field, required = FALSE) {
  value <- row[[field]]
  if (is.null(value) || !length(value) || is.na(value) || value %in% c("", "NA", ".")) {
    if (required) stop("Missing required contrast input: ", field)
    return(NA_real_)
  }
  value <- suppressWarnings(as.numeric(value))
  if (length(value) != 1L || !is.finite(value)) stop("Invalid numeric contrast input: ", field)
  value
}

contrast_rows <- function(leads, comparisons, frequencies, summary) {
  groups <- names(summary$group_counts); totals <- unlist(summary$group_counts)
  if (length(groups) != 3L || anyDuplicated(groups) || any(!is.finite(totals) | totals <= 0 | totals != floor(totals)))
    stop("Exactly three distinct groups with positive integer counts are required.")
  find_group <- function(label) {
    matches <- groups[display_group(groups) == label]
    if (length(matches) != 1L) stop("Expected exactly one source group displayed as ", label, "; groups are never pooled.")
    matches
  }
  control <- find_group("Control")
  cases <- setNames(vapply(c("Diabetes", "Prediabetes"), find_group, character(1)), c("Diabetes", "Prediabetes"))
  if (!all(c("locus_id", "chrom", "lead_pos", "status") %in% names(leads)) || !nrow(leads) || anyDuplicated(leads$locus_id))
    stop("A nonempty, unique selected-locus table is required.")
  observed_status <- table(leads$status); expected_status <- unlist(summary$locus_status_counts)
  if (nrow(leads) != summary$selected_loci || is.null(expected_status) || !setequal(names(observed_status), names(expected_status)) ||
      any(as.numeric(observed_status[names(expected_status)]) != expected_status)) stop("Locus statuses do not match the completed-run summary.")
  required_pairs <- c("locus_id", "haplotype_id", "group_a", "group_b", "test_status", "n_resolved_a", "n_resolved_b", "carriers_a", "carriers_b", "carrier_frequency_a", "carrier_frequency_b")
  if (!all(required_pairs %in% names(comparisons)) || !all(c("locus_id", "haplotype_id", "group", "haplotype_frequency") %in% names(frequencies))) stop("Incomplete contrast input schema.")
  for (data in list(comparisons, frequencies)) {
    if (any(!data$locus_id %in% leads$locus_id)) stop("Input refers to an unknown locus.")
    if (any(leads$status[match(data$locus_id, leads$locus_id)] != "analysed")) stop("An excluded locus cannot contain haplotype results.")
    if (any(is.na(data$haplotype_id) | data$haplotype_id %in% c("", "NA", "."))) stop("Missing haplotype ID.")
  }
  if (any(!comparisons$group_a %in% groups | !comparisons$group_b %in% groups) || any(comparisons$group_a == comparisons$group_b)) stop("Invalid pairwise groups.")
  pair_keys <- paste(comparisons$locus_id, comparisons$haplotype_id, pmin(comparisons$group_a, comparisons$group_b), pmax(comparisons$group_a, comparisons$group_b), sep = "\r")
  if (anyDuplicated(pair_keys)) stop("Duplicate pairwise comparison, including reversed duplicates.")
  if (any(!frequencies$group %in% groups)) stop("Frequency table refers to an unknown group.")
  if (anyDuplicated(paste(frequencies$locus_id, frequencies$haplotype_id, frequencies$group, sep = "\r"))) stop("Duplicate haplotype frequency.")
  for (i in seq_len(nrow(frequencies))) {
    value <- contrast_numeric(as.list(frequencies[i, ]), "haplotype_frequency")
    if (!is.na(value) && (value < 0 || value > 1)) stop("Haplotype frequencies must be within [0,1].")
  }
  output <- list()
  for (i in seq_len(nrow(leads))) {
    lead <- as.list(leads[i, ])
    pairs <- comparisons[comparisons$locus_id == lead$locus_id, , drop = FALSE]
    freq <- frequencies[frequencies$locus_id == lead$locus_id, , drop = FALSE]
    haps <- sort(unique(c(freq$haplotype_id, pairs$haplotype_id)))
    if (!length(haps)) haps <- NA_character_
    for (hap in haps) for (label in names(cases)) {
      case <- cases[[label]]
      row <- setNames(rep(list(NA), length(CONTRAST_FIELDS)), CONTRAST_FIELDS)
      row$locus_id <- lead$locus_id; row$chrom <- lead$chrom; row$lead_pos <- lead$lead_pos
      row$haplotype_id <- hap; row$contrast <- paste(label, "vs Control")
      row$case_group <- case; row$control_group <- control; row$panel_status <- lead$status
      row$comparison_status <- "not_estimable"; row$inference_status <- "not_tested"
      row$reason <- if (lead$status != "analysed") lead$status else if (is.na(hap)) "no_resolved_haplotypes" else "pairwise_result_unavailable"
      if (!is.na(hap) && lead$status == "analysed") {
        source <- pairs[pairs$haplotype_id == hap & ((pairs$group_a == case & pairs$group_b == control) | (pairs$group_a == control & pairs$group_b == case)), , drop = FALSE]
        for (side in c("case", "control")) {
          group <- if (side == "case") case else control
          source_freq <- freq[freq$haplotype_id == hap & freq$group == group, , drop = FALSE]
          if (nrow(source_freq)) row[[paste0("haplotype_frequency_", side, "_descriptive")]] <- contrast_numeric(as.list(source_freq[1, ]), "haplotype_frequency")
        }
        if (nrow(source)) {
          source <- as.list(source[1, ]); forward <- source$group_a == case
          row$reason <- source$test_status; row$inference_status <- source$test_status
          for (side in c("case", "control")) {
            suffix <- if ((side == "case") == forward) "a" else "b"
            for (field in c("n_resolved", "carriers", "carrier_frequency")) row[[paste0(field, "_", side)]] <- contrast_numeric(source, paste0(field, "_", suffix))
            n <- row[[paste0("n_resolved_", side)]]; k <- row[[paste0("carriers_", side)]]; f <- row[[paste0("carrier_frequency_", side)]]
            total <- totals[[if (side == "case") case else control]]
            if (any(!is.na(c(n, k, f))) && anyNA(c(n, k))) stop("Partial participant counts in contrast.")
            if ((!is.na(n) && (n < 0 || n != floor(n) || n > total)) || (!is.na(k) && (k < 0 || k != floor(k) || (!is.na(n) && k > n)))) stop("Invalid participant counts in contrast.")
            if (!is.na(f) && (f < 0 || f > 1)) stop("Invalid carrier frequency.")
            if (!is.na(n) && !is.na(k) && n > 0) {
              expected <- k/n
              if (is.na(f) || abs(f-expected) > 1e-7) stop("Carrier frequency disagrees with participant counts.")
              row[[paste0("carrier_frequency_", side)]] <- expected
            } else if (!is.na(f)) stop("Carrier frequency lacks a positive observed denominator.")
          }
          observed <- vapply(c("n_resolved_case", "n_resolved_control", "carriers_case", "carriers_control"), function(field) row[[field]], numeric(1))
          if (!anyNA(observed) && all(observed[1:2] > 0)) {
            row$descriptive_carrier_frequency_difference <- row$carrier_frequency_case-row$carrier_frequency_control
            row$comparison_status <- "descriptive_only"
          }
          if (source$test_status == "tested_unadjusted") {
            params <- summary$parameters
            if (!isTRUE(params$test_associations) || !isTRUE(params$reference_genotypes_confirmed)) stop("Tested effects require verified reference-call provenance.")
            if (row$comparison_status != "descriptive_only") stop("Tested contrast lacks valid participant counts.")
            if (is.null(params$min_test_samples) || is.null(params$min_resolved_fraction) || is.null(params$min_carriers)) stop("Missing association eligibility thresholds.")
            for (side in c("case", "control")) {
              n <- row[[paste0("n_resolved_", side)]]; total <- totals[[if (side == "case") case else control]]
              if (n < params$min_test_samples || n/total < params$min_resolved_fraction) stop("Tested contrast fails participant-resolution requirements.")
            }
            n <- sum(observed[1:2]); k <- sum(observed[3:4])
            if (min(k, n-k) < params$min_carriers) stop("Tested contrast fails carrier-count requirements.")
            difference <- contrast_numeric(source, "carrier_frequency_difference_a_minus_b", TRUE)
            low <- contrast_numeric(source, "carrier_difference_CI95_low", TRUE); high <- contrast_numeric(source, "carrier_difference_CI95_high", TRUE)
            p <- contrast_numeric(source, "p_fisher_two_sided", TRUE); q <- contrast_numeric(source, "q_BH_all_pairwise_tests", TRUE)
            if (!forward) { difference <- -difference; old_low <- low; low <- -high; high <- -old_low }
            if (low > difference || high < difference || low < -1-1e-8 || high > 1+1e-8 || abs(difference-row$descriptive_carrier_frequency_difference) > 1e-7 || p < 0 || q < 0 || q > 1 || p > q+1e-8) stop("Invalid tested effect, interval or p/q values.")
            row$comparison_status <- "tested"; row$carrier_frequency_difference <- difference
            row$CI95_low <- low; row$CI95_high <- high; row$p_fisher_two_sided <- p; row$q_BH_all_pairwise_tests <- q
          } else {
            for (field in CONTRAST_INFERENCE_FIELDS)
              if (!is.na(contrast_numeric(source, field))) stop("Untested source contains inferential values: ", field)
          }
        }
      }
      output[[length(output)+1L]] <- row
    }
  }
  output
}

contrast_reason_label <- function(reason, cutoff) {
  labels <- c(lead_below_min_call_rate = sprintf("Lead genotypes below %.0f%% completeness", cutoff*100),
    palindromic_lead_orientation_unconfirmed = "Palindromic lead orientation unconfirmed",
    negative_liftover_allele_orientation_unconfirmed = "Lead-allele orientation unconfirmed",
    no_unique_PASS_sequence_allele_match_at_lead = "No unique matched lead record", multiple_matching_lead_records = "Multiple matching lead records",
    fewer_than_two_eligible_markers = "Fewer than two eligible markers", not_tested_reference_provenance = "Reference-call provenance unconfirmed",
    too_few_resolved_samples = "Too few resolved participants", low_resolved_fraction = "Insufficient phase resolution",
    no_carrier_variation = "No carrier variation for the test", too_few_carriers_or_noncarriers = "Too few carriers or noncarriers",
    not_requested = "Association test not requested", pairwise_result_unavailable = "Pairwise result unavailable", no_resolved_haplotypes = "No resolved haplotypes")
  if (reason %in% names(labels)) labels[[reason]] else gsub("_", " ", reason, fixed = TRUE)
}

contrast_device <- function(path, type, width, height, dpi = 220) {
  if (type == "png") png(path, width = width, height = height, units = "in", res = dpi, type = "cairo")
  else if (type == "pdf") pdf(path, width = width, height = height, useDingbats = FALSE)
  else if (type == "svg") svg(path, width = width, height = height)
  else stop("Unsupported figure format.")
  par(mar = rep(0, 4), family = "sans", bg = "white")
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
}

contrast_interpretation <- function(rows, summary) {
  # Region status placeholders have no identified haplotype and are not
  # unavailable comparisons for actual haplotypes.
  has_haplotype <- function(r) !is.na(r$haplotype_id) && !r$haplotype_id %in% c("", "NA", ".")
  haplotype_rows <- Filter(has_haplotype, rows)
  placeholders <- Filter(function(r) !has_haplotype(r), rows)
  excluded <- Filter(function(r) r$panel_status != "analysed", placeholders)
  unresolved <- Filter(function(r) r$panel_status == "analysed", placeholders)
  tested <- Filter(function(r) r$comparison_status == "tested", haplotype_rows)
  descriptive <- Filter(function(r) r$comparison_status == "descriptive_only", haplotype_rows)
  unavailable <- Filter(function(r) r$comparison_status == "not_estimable", haplotype_rows)
  observed <- Filter(function(r) !is.na(r$descriptive_carrier_frequency_difference), haplotype_rows)
  region_count <- function(data) length(unique(vapply(data, `[[`, character(1), "locus_id")))
  patterns <- length(unique(vapply(observed, function(r) paste(r$locus_id, r$haplotype_id, sep = "\r"), character(1))))
  eligible <- region_count(Filter(function(r) r$panel_status == "analysed", rows))
  source_eligible <- summary$locus_status_counts$analysed
  if (is.null(source_eligible)) source_eligible <- 0L
  lines <- c("Haplotype participant-carrier contrasts against Control", "",
    "The source analysis completed. Unavailable inference reflects analysis gates, not a plotting failure.",
    "Requested contrasts: Diabetes vs Control and Prediabetes vs Control.",
    sprintf("Report scope: %d of %d selected regions; %d analysed regions in this report. The complete source has %d analysed (panel-eligible) regions.", region_count(rows), summary$selected_loci, eligible, source_eligible),
    sprintf("Observed carrier differences are available for %d local haplotypes in %d regions.", patterns, region_count(observed)),
    sprintf("%d tested haplotype comparisons; %d descriptive-only differences; %d identified-haplotype comparisons lack a computable carrier contrast.", length(tested), length(descriptive), length(unavailable)),
    sprintf("%d excluded-region placeholder rows cover %d excluded regions in this report. These rows record requested contrasts without an identified haplotype; they are not haplotype comparisons.", length(excluded), region_count(excluded)),
    "Observed differences equal carriers/resolved participants in Diabetes or Prediabetes minus the corresponding Control fraction.",
    "These are participant-carrier differences, distinct from chromosome-copy haplotype frequencies.",
    "Observed descriptive differences condition on callable, phase-resolved people; missingness and unconfirmed reference-call provenance can bias them.",
    "Descriptive-only differences have no confidence interval, p-value, q-value or significance claim.",
    "Tested rows alone retain the source Newcombe-Wilson 95% intervals and full original pairwise BH correction; this plot performs no inference.")
  if (length(unresolved)) lines <- c(lines, sprintf("%d further placeholder rows cover %d analysed regions with no resolved haplotypes; these are not identified-haplotype comparisons.", length(unresolved), region_count(unresolved)))
  if (identical(summary$parameters$reference_genotypes_confirmed, FALSE)) lines <- c(lines,
    "Source parameter reference_genotypes_confirmed=false: reference-call provenance is unconfirmed, so the reference-provenance gate blocks association inference.",
    "Confirming reference provenance alone does not establish eligibility; participant-resolution and carrier-count requirements must also pass. The reported reason is the first failing gate; other limitations may also apply.")
  reference_calls <- summary$scan_counts$high_quality_homozygous_reference_genotypes
  if (is.numeric(reference_calls) && length(reference_calls) == 1L && !is.na(reference_calls) && is.finite(reference_calls) && reference_calls >= 0 && reference_calls == floor(reference_calls)) lines <- c(lines,
    sprintf("Source scan recorded %.0f high-quality homozygous-reference genotypes in the selected windows. This post-QC count does not establish whether reference calls were absent before QC.", reference_calls))
  if (!length(tested)) lines <- c(lines, "No inferential haplotype association can be established from this run. Available observed differences are shown without significance claims.")
  for (label in c("Diabetes vs Control", "Prediabetes vs Control")) {
    selected <- Filter(function(r) r$contrast == label, observed)
    if (!length(selected)) next
    delta <- vapply(selected, `[[`, numeric(1), "descriptive_carrier_frequency_difference")
    best <- selected[[which.max(abs(delta))]]
    lines <- c(lines, sprintf("Largest absolute observed %s difference: %s, %+.2f percentage points (%d/%d vs %d/%d carriers/resolved people); descriptive ranking, not evidence of significance.", label, best$haplotype_id, 100*best$descriptive_carrier_frequency_difference, best$carriers_case, best$n_resolved_case, best$carriers_control, best$n_resolved_control))
  }
  if (length(descriptive)) {
    reasons <- table(vapply(descriptive, `[[`, character(1), "inference_status"))
    lines <- c(lines, "", "Reasons observed haplotype comparisons remain descriptive-only (comparison counts; first failing gate; other limitations may also apply):", paste(as.integer(reasons), names(reasons), sep = ": "))
  }
  if (length(unavailable)) {
    reasons <- table(vapply(unavailable, `[[`, character(1), "reason"))
    lines <- c(lines, "", "Reasons identified-haplotype carrier contrasts are unavailable (comparison counts):", paste(as.integer(reasons), names(reasons), sep = ": "))
  }
  if (length(excluded)) {
    by_region <- excluded[!duplicated(vapply(excluded, `[[`, character(1), "locus_id"))]
    reasons <- table(vapply(by_region, `[[`, character(1), "reason"))
    lines <- c(lines, "", "Excluded-region reasons (region counts, not haplotype comparison counts):", paste(as.integer(reasons), names(reasons), sep = ": "))
  }
  c(lines, "", "Unavailable means missing, not zero difference or equal groups.",
    "Causal variants remain unresolved; observed contrasts and external GWAS associations do not establish causality.",
    "Recover validated reference/alternate genotypes and usable phase before attempting gated inference. All results are private artifacts.")
}

# Compact overview: all observed haplotypes, one facet per region and two columns.
draw_haplotype_contrasts <- function(rows, summary, path, type = "png", columns = 3L) {
  observed <- Filter(function(r) !is.na(r$descriptive_carrier_frequency_difference), rows)
  loci <- unique(vapply(observed, `[[`, character(1), "locus_id"))
  tested <- sum(vapply(rows, function(r) r$comparison_status == "tested", logical(1)))
  descriptive <- sum(vapply(rows, function(r) r$comparison_status == "descriptive_only", logical(1)))
  if (!length(loci)) {
    # Preserve status-only rendering when no observed carrier difference exists.
    n <- length(rows)/2; height <- max(7, 3.5+.30*n)
    contrast_device(path, type, 14, height); on.exit(dev.off())
    text(.04, .96, "Haplotype contrasts against Control", adj = 0, font = 2, cex = 1.6, col = "#172f46")
    text(.04, .915, sprintf("%d selected regions | No observed carrier differences available", summary$selected_loci), adj = 0, col = "#526579")
    text(c(.49, .82), .85, c("Diabetes vs Control", "Prediabetes vs Control"), font = 2, col = c("#b64b18", "#008666"))
    y <- seq(.79, .20, length.out = n)
    for (i in seq_len(n)) {
      pair <- rows[(2*i-1):(2*i)]; text(.29, y[i], pair[[1]]$locus_id, adj = 1, cex = .8)
      for (j in 1:2) {
        x <- c(.49, .82)[j]; rect(x-.15, y[i]-.25/max(n, 1), x+.15, y[i]+.25/max(n, 1), col = "#f2f5f8", border = "white")
        text(x, y[i], paste("Not estimable", contrast_reason_label(pair[[j]]$reason, summary$parameters$min_call_rate), sep = "\n"), cex = .65, col = "#526579")
      }
    }
    text(.04, .11, "No observed differences can be shown; unavailable comparisons do not mean zero effects or equal groups.", adj = 0, cex = .88)
    text(.04, .06, "Private outputs | Independent R | Causal variants unresolved", adj = 0, cex = .8, col = "#63758a")
    return(invisible(NULL))
  }
  by_row <- split(seq_along(loci), ceiling(seq_along(loci)/columns))
  max_haps <- vapply(by_row, function(ix) max(vapply(loci[ix], function(locus) length(unique(vapply(Filter(function(r) r$locus_id == locus, rows), `[[`, character(1), "haplotype_id"))), integer(1))), integer(1))
  row_heights <- pmax(2.4, 1.05+.24*max_haps); height <- 2.5+sum(row_heights); width <- 4.6*columns
  contrast_device(path, type, width, height); on.exit(dev.off())
  text(.025, 1-.24/height, "Observed haplotype carrier differences against Control", adj = 0, font = 2, cex = 1.5, col = "#172f46")
  text(.025, 1-.55/height, sprintf("%d observed regions | %d tested + %d descriptive-only differences | %d selected regions without observed differences", length(loci), tested, descriptive, summary$selected_loci-length(loci)), adj = 0, cex = .86, col = "#526579")
  text(.025, 1-.83/height, "Signed percentage points: carriers / resolved people in disease group minus Control. Descriptive cells have no CI or significance test.", adj = 0, cex = .84, col = "#526579")
  scale <- max(abs(vapply(observed, `[[`, numeric(1), "descriptive_carrier_frequency_difference"))); if (scale == 0) scale <- .01
  palette <- colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(201)
  gap <- .032; margin <- .025; facet_width <- (1-2*margin-gap*(columns-1))/columns
  cumulative <- c(0, cumsum(row_heights))/height; top <- 1-1.20/height
  for (i in seq_along(loci)) {
    locus <- loci[i]; grid_row <- (i-1L)%/%columns+1L; grid_col <- (i-1L)%%columns+1L
    left <- margin+(grid_col-1)*(facet_width+gap); right <- left+facet_width
    facet_top <- top-cumulative[grid_row]; facet_bottom <- top-cumulative[grid_row+1L]
    local <- Filter(function(r) r$locus_id == locus, rows); haps <- unique(vapply(local, `[[`, character(1), "haplotype_id")); first <- local[[1]]
    text(left, facet_top-.1/height, paste0(locus, "  ", first$chrom, ":", first$lead_pos), adj = 0, font = 2, cex = .91, col = "#172f46")
    text(left, facet_top-.34/height, "Observed carrier difference (percentage points)", adj = 0, cex = .64, col = "#526579")
    matrix_left <- left+.065; cell_width <- (right-matrix_left)/2; x <- matrix_left+(1:2-.5)*cell_width
    text(x, facet_top-.58/height, c("Diabetes - Control", "Prediabetes - Control"), cex = .61, col = "#40566d")
    matrix_top <- facet_top-.76/height; matrix_bottom <- facet_bottom+.16/height; step <- (matrix_top-matrix_bottom)/length(haps)
    for (j in seq_along(haps)) {
      hap <- haps[j]; y <- matrix_top-(j-.5)*step
      text(matrix_left-.007, y, sub(paste0(locus, "_"), "", hap, fixed = TRUE), adj = 1, cex = .69, col = "#233445")
      for (k in 1:2) {
        cell <- Filter(function(r) r$haplotype_id == hap && r$contrast == c("Diabetes vs Control", "Prediabetes vs Control")[k], local)[[1]]
        value <- cell$descriptive_carrier_frequency_difference
        color <- if (is.na(value)) "#d9dde3" else palette[round((value/scale+1)*100)+1L]
        rect(x[k]-.48*cell_width, y-.46*step, x[k]+.48*cell_width, y+.46*step, col = color, border = NA)
        label <- if (is.na(value)) "NA" else sprintf("%+.2f%s", value*100, if (cell$comparison_status == "tested") "*" else "")
        text(x[k], y, label, cex = .77, col = if (!is.na(value) && abs(value)/scale > .65) "white" else "#233445")
      }
    }
  }
  ly <- .69/height
  for (k in seq_len(201)) rect(.07+(k-1)*.0015, ly, .07+k*.0015, ly+.15/height, col = palette[k], border = NA)
  text(c(.07, .22, .3715), ly-.12/height, c(sprintf("%.1f pp", -scale*100), "0", sprintf("%+.1f pp", scale*100)), cex = .72, col = "#526579")
  text(.42, ly+.07/height, "* Source-tested row; see regional CI and q. Gray = unavailable; full exclusions in contrast_data.tsv.", adj = 0, cex = .67, col = "#526579")
  text(.025, .16/height, paste(if (isTRUE(summary$synthetic)) "SYNTHETIC" else "Real cohort", "| Independent R | Private results | Selected observed calls can bias differences | Causal variants unresolved"), adj = 0, cex = .73, col = "#63758a")
}

# One native R forest-style figure per observed region. Descriptive points have no CI.
draw_regional_haplotype_contrasts <- function(rows, summary, path, type = "png") {
  loci <- unique(vapply(rows, `[[`, character(1), "locus_id")); if (length(loci) != 1L) stop("Regional contrasts require one locus.")
  haps <- unique(vapply(rows, `[[`, character(1), "haplotype_id")); height <- max(5.5, 2.6+.47*length(haps))
  contrast_device(path, type, 14, height); on.exit(dev.off())
  tested <- any(vapply(rows, function(r) r$comparison_status == "tested", logical(1)))
  text(.035, .947, paste(loci, "haplotype carrier contrasts"), adj = 0, font = 2, cex = 1.5, col = "#172f46")
  text(.035, .889, "Observed carrier fractions among resolved people; positive differences indicate more carriers than in controls.", adj = 0, cex = .85, col = "#526579")
  regions <- list(c(.14, .49), c(.60, .95)); colors <- c("#b64b18", "#008666")
  contrasts <- c("Diabetes vs Control", "Prediabetes vs Control"); top <- .73; bottom <- .26; step <- (top-bottom)/length(haps)
  bounds <- unlist(lapply(rows, function(r) c(r$descriptive_carrier_frequency_difference, if (r$comparison_status == "tested") c(r$CI95_low, r$CI95_high) else NULL)))
  limit <- max(.05, abs(bounds[is.finite(bounds)]))*1.14
  for (k in 1:2) {
    region <- regions[[k]]; plot_left <- region[1]+.005; plot_right <- region[2]-.14
    position <- function(d) plot_left+(d+limit)/(2*limit)*(plot_right-plot_left)
    text(mean(region), .81, contrasts[k], font = 2, col = colors[k], cex = 1.07)
    segments(position(0), bottom, position(0), top, lty = 3, col = "#8090a0")
    for (j in seq_along(haps)) {
      y <- top-(j-.5)*step; cell <- Filter(function(r) r$haplotype_id == haps[j] && r$contrast == contrasts[k], rows)[[1]]
      rect(region[1], y-.46*step, region[2], y+.46*step, col = if (j%%2) "#f3f6f8" else "white", border = NA)
      segments(position(0), y-.45*step, position(0), y+.45*step, lty = 3, col = "#8090a0")
      if (k == 1) text(region[1]-.014, y, sub(paste0(loci, "_"), "", haps[j], fixed = TRUE), adj = 1, cex = .8)
      effect <- cell$descriptive_carrier_frequency_difference
      if (is.na(effect)) text(mean(region), y, "Not estimable", cex = .73, col = "#526579") else {
        if (cell$comparison_status == "tested") segments(position(cell$CI95_low), y, position(cell$CI95_high), y, col = colors[k], lwd = 1.6)
        points(position(effect), y, pch = if (cell$comparison_status == "tested") 16 else 1, cex = .86, col = colors[k])
        detail <- sprintf("%+.2f pp | %d/%d vs %d/%d", 100*effect, cell$carriers_case, cell$n_resolved_case, cell$carriers_control, cell$n_resolved_control)
        text(region[2]-.003, y+.1*step, detail, adj = 1, cex = .63, col = "#40566d")
        text(region[2]-.003, y-.19*step, if (cell$comparison_status == "tested") sprintf("95%% CI; BH q=%.3g", cell$q_BH_all_pairwise_tests) else "Descriptive only; inference withheld", adj = 1, cex = .54, col = "#63758a")
      }
    }
    ticks <- pretty(c(-limit, limit), n = 4); ticks <- ticks[abs(ticks) <= limit]
    text(vapply(ticks, position, numeric(1)), bottom-.025, sprintf("%+.0f", 100*ticks), cex = .66, col = "#526579")
    text(mean(c(plot_left, plot_right)), bottom-.066, "Carrier difference (percentage points)", cex = .69, col = "#526579")
  }
  text(.035, .114, if (tested) "Filled points: tested with source 95% CI. Open points: observed differences only, without inferential uncertainty." else "Open points show descriptive differences only. No confidence intervals, p-values or significance claims are available.", adj = 0, cex = .81, col = "#40566d")
  text(.035, .070, "Missingness and unconfirmed reference-call provenance can bias these observed fractions. Causal variants remain unresolved.", adj = 0, cex = .78, col = "#526579")
  text(.035, .031, "Private project results | Independent native R | Denominators: resolved participants", adj = 0, cex = .72, col = "#63758a")
}

contrast_sha256 <- function(paths) {
  command <- Sys.which("sha256sum"); if (!nzchar(command)) stop("sha256sum is required for contrast input provenance.")
  values <- vapply(paths, function(path) {
    value <- system2(command, shQuote(normalizePath(path)), stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(value, "status")) || length(value) != 1L || !grepl("^[a-f0-9]{64} ", value)) stop("Could not fingerprint contrast input.")
    substr(value, 1, 64)
  }, character(1))
  setNames(unname(values), basename(paths))
}

export_haplotype_contrasts <- function(run_dir, out) {
  if (!file.exists(file.path(run_dir, "COMPLETED.txt"))) stop("A completed analysis is required.")
  if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) stop("Use a new or empty output directory.")
  paths <- file.path(run_dir, CONTRAST_SOURCES); before <- contrast_sha256(paths)
  summary <- jsonlite::fromJSON(paths[1], simplifyVector = FALSE)
  if (!isTRUE(summary$inputs_stable_during_read)) stop("Source analysis did not verify stable inputs.")
  read_input <- function(path) {
    header <- strsplit(readLines(path, n = 1L), "\t", fixed = TRUE)[[1]]
    if (!length(header) || anyDuplicated(header)) stop("Missing or duplicate table headers.")
    read_tsv(path)
  }
  rows <- contrast_rows(read_input(paths[2]), read_input(paths[3]), read_input(paths[4]), summary)
  if (!identical(before, contrast_sha256(paths))) stop("Source analysis changed while reading contrast inputs.")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write_rows(file.path(out, "contrast_data.tsv"), rows, CONTRAST_FIELDS)
  for (type in c("png", "pdf", "svg")) draw_haplotype_contrasts(rows, summary, file.path(out, paste0("haplotype_contrasts.", type)), type)
  writeLines(contrast_interpretation(rows, summary), file.path(out, "INTERPRETATION.txt"))
  observed <- Filter(function(r) !is.na(r$descriptive_carrier_frequency_difference), rows)
  loci <- unique(vapply(observed, `[[`, character(1), "locus_id"))
  for (locus in loci) {
    regional <- Filter(function(r) r$locus_id == locus, rows); folder <- file.path(out, "regional", locus)
    if (!grepl("^[A-Za-z0-9_-]+$", locus)) stop("Unsafe locus identifier for an output path.")
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    for (type in c("png", "pdf")) draw_regional_haplotype_contrasts(regional, summary, file.path(folder, paste0("haplotype_contrasts.", type)), type)
    writeLines(contrast_interpretation(regional, summary), file.path(folder, "INTERPRETATION.txt"))
  }
  index <- c("# Haplotype contrasts: every selected region", "",
    sprintf("Source completed analysis: `%s`.", normalizePath(run_dir)), "",
    "Observed participant-carrier differences are descriptive unless a row is explicitly tested. All absent estimates retain their exclusion reason.", "",
    "[Overview](haplotype_contrasts.png) | [Full aggregate table](contrast_data.tsv) | [Interpretation](INTERPRETATION.txt)", "",
    "| Region | Panel status | Observed disease-control contrasts | Regional picture | Interpretation |",
    "| --- | --- | ---: | --- | --- |")
  for (locus in unique(vapply(rows, `[[`, character(1), "locus_id"))) {
    local <- Filter(function(r) r$locus_id == locus, rows)
    n_observed <- sum(vapply(local, function(r) !is.na(r$descriptive_carrier_frequency_difference), logical(1)))
    index <- c(index, sprintf("| %s | %s | %d | %s | %s |", locus, local[[1]]$panel_status, n_observed,
      if (n_observed) sprintf("[PNG](regional/%s/haplotype_contrasts.png), [PDF](regional/%s/haplotype_contrasts.pdf)", locus, locus) else "Unavailable",
      if (n_observed) sprintf("[Read](regional/%s/INTERPRETATION.txt)", locus) else local[[1]]$reason))
  }
  writeLines(index, file.path(out, "index.md"))
  if (!identical(before, contrast_sha256(paths))) stop("Source analysis changed during contrast plotting; discard these outputs.")
  statuses <- table(factor(vapply(rows, `[[`, character(1), "comparison_status"), levels = c("tested", "descriptive_only", "not_estimable")))
  manifest <- list(private_results = TRUE, implementation = "independent_R", source_run = basename(normalizePath(run_dir)),
    source_run_path = normalizePath(run_dir),
    source_implementation = summary$implementation %||% "Python", source_sha256 = as.list(before),
    inputs_stable_during_plotting = TRUE, synthetic = isTRUE(summary$synthetic), R_version = R.version.string,
    jsonlite_version = as.character(packageVersion("jsonlite")), rows = length(rows), selected_regions = summary$selected_loci,
    observed_regions = length(loci), comparison_status_counts = as.list(setNames(as.integer(statuses), names(statuses))),
    created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), png_dpi = 220,
    inference_performed_by_renderer = FALSE, quantity = "participant-carrier frequency difference; descriptive values separate from tested inference")
  jsonlite::write_json(manifest, file.path(out, "manifest.json"), pretty = TRUE, auto_unbox = TRUE)
  cat("Native R carrier contrasts: ", length(rows), " rows; ", length(loci), " observed regions; ", statuses[["tested"]], " tested; ", statuses[["descriptive_only"]], " descriptive only.\n", sep = "")
  invisible(rows)
}

if (sys.nframe() == 0L) {
  argv <- commandArgs(TRUE)
  if (length(argv) != 4L || !setequal(argv[c(1, 3)], c("--run-dir", "--out"))) stop("Use: Rscript code/R/plot_haplotype_contrasts.R --run-dir COMPLETED_RUN --out NEW_DIRECTORY")
  args <- as.list(argv[c(2, 4)]); names(args) <- gsub("-", "_", sub("^--", "", argv[c(1, 3)]))
  do.call(export_haplotype_contrasts, args)
}
