#!/usr/bin/env Rscript
# Native R disease-versus-control contrast figure from completed aggregate tables.
# Missing/unavailable estimates remain NA and are never drawn at zero.
entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
if (!exists("read_tsv", mode = "function")) source(file.path(dirname(normalizePath(entry)), "compare_haplotypes.R"))

CONTRAST_FIELDS <- strsplit("locus_id chrom lead_pos haplotype_id contrast case_group control_group panel_status comparison_status reason n_resolved_case n_resolved_control carriers_case carriers_control carrier_frequency_case carrier_frequency_control carrier_frequency_difference CI95_low CI95_high p_fisher_two_sided q_BH_all_pairwise_tests haplotype_frequency_case_descriptive haplotype_frequency_control_descriptive", " ")[[1]]

contrast_rows <- function(leads, comparisons, frequencies, summary) {
  groups <- names(summary$group_counts)
  find_group <- function(label) {
    matches <- groups[display_group(groups) == label]
    if (length(matches) != 1L) stop("Expected exactly one source group displayed as ", label, "; groups are never pooled.")
    matches
  }
  control <- find_group("Control")
  cases <- setNames(vapply(c("Diabetes", "Prediabetes"), find_group, character(1)), c("Diabetes", "Prediabetes"))
  numeric_value <- function(row, column) {
    value <- row[[column]]
    if (is.null(value) || !length(value) || value %in% c("", "NA")) return(NA_real_)
    value <- suppressWarnings(as.numeric(value))
    if (!is.finite(value)) stop("Invalid numeric contrast input: ", column)
    value
  }
  output <- list()
  for (i in seq_len(nrow(leads))) {
    lead <- as.list(leads[i, ])
    pairs <- comparisons[comparisons$locus_id == lead$locus_id, , drop = FALSE]
    freq <- frequencies[frequencies$locus_id == lead$locus_id, , drop = FALSE]
    haps <- if (lead$status == "analysed") sort(unique(c(freq$haplotype_id, pairs$haplotype_id))) else character()
    haps <- haps[!is.na(haps) & !haps %in% c("", "NA")]
    if (!length(haps)) haps <- NA_character_
    for (hap in haps) for (label in names(cases)) {
      case <- cases[[label]]
      row <- setNames(rep(list(NA), length(CONTRAST_FIELDS)), CONTRAST_FIELDS)
      row$locus_id <- lead$locus_id; row$chrom <- lead$chrom; row$lead_pos <- lead$lead_pos
      row$haplotype_id <- hap; row$contrast <- paste(label, "vs Control")
      row$case_group <- case; row$control_group <- control; row$panel_status <- lead$status
      row$comparison_status <- "not_estimable"
      row$reason <- if (lead$status != "analysed") lead$status else if (is.na(hap)) "no_resolved_haplotypes" else "pairwise_result_unavailable"
      if (!is.na(hap) && lead$status == "analysed") {
        source <- pairs[pairs$haplotype_id == hap & ((pairs$group_a == case & pairs$group_b == control) |
                        (pairs$group_a == control & pairs$group_b == case)), , drop = FALSE]
        if (nrow(source) > 1L) stop("Duplicate disease/control pairwise result.")
        for (side in c("case", "control")) {
          group <- if (side == "case") case else control
          source_freq <- freq[freq$haplotype_id == hap & freq$group == group, , drop = FALSE]
          if (nrow(source_freq) > 1L) stop("Duplicate haplotype/group frequency.")
          if (nrow(source_freq)) row[[paste0("haplotype_frequency_", side, "_descriptive")]] <- numeric_value(as.list(source_freq[1, ]), "haplotype_frequency")
        }
        if (nrow(source)) {
          source <- as.list(source[1, ])
          forward <- source$group_a == case
          row$reason <- source$test_status
          for (side in c("case", "control")) {
            suffix <- if ((side == "case") == forward) "a" else "b"
            for (field in c("n_resolved", "carriers", "carrier_frequency")) row[[paste0(field, "_", side)]] <- numeric_value(source, paste0(field, "_", suffix))
          }
          if (source$test_status == "tested_unadjusted") {
            difference <- numeric_value(source, "carrier_frequency_difference_a_minus_b")
            low <- numeric_value(source, "carrier_difference_CI95_low")
            high <- numeric_value(source, "carrier_difference_CI95_high")
            p <- numeric_value(source, "p_fisher_two_sided")
            q <- numeric_value(source, "q_BH_all_pairwise_tests")
            if (anyNA(c(difference, low, high, p, q)) || low > difference || high < difference ||
                any(c(difference, low, high) < -1 | c(difference, low, high) > 1) || any(c(p, q) < 0 | c(p, q) > 1))
              stop("Tested contrast lacks a valid effect, interval, p-value or BH q-value.")
            row$comparison_status <- "tested"; row$reason <- source$test_status
            row$carrier_frequency_difference <- if (forward) difference else -difference
            row$CI95_low <- if (forward) low else -high
            row$CI95_high <- if (forward) high else -low
            row$p_fisher_two_sided <- p; row$q_BH_all_pairwise_tests <- q
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
    no_unique_PASS_sequence_allele_match_at_lead = "No unique matched lead record",
    multiple_matching_lead_records = "Multiple matching lead records",
    fewer_than_two_eligible_markers = "Fewer than two eligible markers",
    not_tested_reference_provenance = "Reference-call provenance unconfirmed",
    too_few_resolved_samples = "Too few resolved participants", low_resolved_fraction = "Insufficient phase resolution",
    no_carrier_variation = "No carrier variation for the test",
    too_few_carriers_or_noncarriers = "Too few carriers or noncarriers",
    not_requested = "Association test not requested", pairwise_result_unavailable = "Pairwise result unavailable",
    no_resolved_haplotypes = "No resolved haplotypes")
  if (reason %in% names(labels)) labels[[reason]] else gsub("_", " ", reason, fixed = TRUE)
}

draw_haplotype_contrasts <- function(rows, summary, path, type = "png") {
  key <- vapply(rows, function(r) paste(r$locus_id, if (is.na(r$haplotype_id)) "" else r$haplotype_id, sep = ":"), character(1))
  keys <- unique(key); n <- length(keys)
  if (!n) stop("No selected loci are available for the contrast figure.")
  height <- max(7, 3.5 + .38*n)
  if (type == "png") png(path, width = 14, height = height, units = "in", res = 300, type = "cairo")
  else if (type == "pdf") pdf(path, width = 14, height = height, useDingbats = FALSE)
  else if (type == "svg") svg(path, width = 14, height = height)
  else stop("Unsupported figure format.")
  on.exit(dev.off())
  par(mar = rep(0, 4), family = "sans", bg = "white")
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  text(.04, .957, "Diabetes and prediabetes haplotypes vs controls", adj = 0, font = 2, cex = 1.6, col = "#172f46")
  counts <- unlist(summary$group_counts); labels <- display_group(names(counts))
  npeople <- sum(counts); tested <- sum(vapply(rows, function(r) r$comparison_status == "tested", logical(1)))
  subtitle <- sprintf("%d participants | %d selected regions | %d estimable disease-control contrasts", npeople, summary$selected_loci, tested)
  text(.04, .915, subtitle, adj = 0, cex = .95, col = "#526579")
  text(.04, .884, if (tested) "Participant-carrier frequency differences and 95% confidence intervals; q values retain the full pairwise BH family." else
    "No eligible haplotype panels: diabetes-control and prediabetes-control differences cannot be estimated.", adj = 0, cex = .9, col = "#526579")
  regions <- list(c(.31, .62), c(.655, .965)); names(regions) <- c("Diabetes vs Control", "Prediabetes vs Control")
  colors <- c("#d55e00", "#009e73")
  control_n <- counts[labels == "Control"]
  for (j in seq_along(regions)) {
    label <- sub(" vs Control", "", names(regions)[j], fixed = TRUE)
    text(mean(regions[[j]]), .825, names(regions)[j], font = 2, col = colors[j], cex = 1.1)
    text(mean(regions[[j]]), .797, sprintf("n = %d vs %d", counts[labels == label], control_n), cex = .8, col = "#526579")
  }
  top <- .763; bottom <- .195; step <- (top-bottom)/n
  y <- top - (seq_len(n)-.5)*step
  if (tested) for (region in regions) segments(mean(region), bottom, mean(region), top, lty = 3, col = "#aab6c3")
  for (i in seq_along(keys)) {
    entries <- rows[key == keys[i]]; first <- entries[[1]]
    label <- paste(first$locus_id, paste0(first$chrom, ":", format(as.numeric(first$lead_pos), big.mark = ",", scientific = FALSE, trim = TRUE)))
    if (!is.na(first$haplotype_id)) label <- paste(first$locus_id, sub(paste0(first$locus_id, "_"), "", first$haplotype_id, fixed = TRUE))
    text(.285, y[i], label, adj = 1, cex = .77, col = "#233445")
    for (j in seq_along(regions)) {
      cell <- Filter(function(r) r$contrast == names(regions)[j], entries)
      if (length(cell) != 1L) stop("Each row requires both disease-control contrasts.")
      cell <- cell[[1]]; region <- regions[[j]]
      background <- if (cell$reason == "palindromic_lead_orientation_unconfirmed") "#fff0ca" else if (i %% 2L) "#f2f5f8" else "#fafbfd"
      rect(region[1], y[i]-.46*step, region[2], y[i]+.46*step, col = background, border = "white")
      if (cell$comparison_status == "tested") {
        # Reserve the rightmost section of each column for BH q labels.
        left <- region[1]+.025; right <- region[2]-.075
        position <- function(d) left+(d+1)/2*(right-left)
        segments(position(0), y[i]-.45*step, position(0), y[i]+.45*step, lty = 3, col = "#aab6c3")
        segments(position(cell$CI95_low), y[i], position(cell$CI95_high), y[i], col = colors[j], lwd = 1.6)
        points(position(cell$carrier_frequency_difference), y[i], pch = 16, cex = .7, col = colors[j])
        text(region[2]-.01, y[i], sprintf("q = %.3g", cell$q_BH_all_pairwise_tests), adj = 1, cex = .62, col = "#526579")
      } else {
        text(mean(region), y[i]+.14*step, "Not estimable", font = 2, cex = .77, col = "#526579")
        reason <- contrast_reason_label(cell$reason, summary$parameters$min_call_rate)
        text(mean(region), y[i]-.19*step, reason, cex = .57, col = "#63758a")
      }
    }
  }
  if (tested) {
    for (region in regions) {
      left <- region[1]+.025; right <- region[2]-.075
      for (d in c(-1, -.5, 0, .5, 1)) text(left+(d+1)/2*(right-left), bottom-.014, sprintf("%+.0f", d*100), cex = .64, col = "#526579")
      text(mean(c(left, right)), bottom-.041, "Carrier-frequency difference (percentage points)", cex = .62, col = "#526579")
    }
  }
  note <- if (tested) "Positive differences indicate more carriers in diabetes/prediabetes than controls; these are exploratory associations." else
    "No haplotype difference can be identified from these inputs. Unavailable estimates are not zero effects."
  text(.04, .116, note, adj = 0, cex = .86, col = "#40566d")
  text(.04, .080, "Causal variants remain unresolved; haplotype or external GWAS associations alone do not establish causality.", adj = 0, cex = .82, col = "#526579")
  text(.04, .04, paste(if (isTRUE(summary$synthetic)) "SYNTHETIC SOFTWARE FIXTURE" else "Real cohort", "| Independent R | Missing genotypes retained | Private project outputs"), adj = 0, cex = .72, col = "#63758a")
}

export_haplotype_contrasts <- function(run_dir, out) {
  if (!file.exists(file.path(run_dir, "COMPLETED.txt"))) stop("A completed analysis is required.")
  if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) stop("Use a new or empty output directory.")
  summary <- jsonlite::fromJSON(file.path(run_dir, "run_summary.json"), simplifyVector = FALSE)
  leads <- read_tsv(file.path(run_dir, "lead_loci.tsv"))
  comparisons <- read_tsv(file.path(run_dir, "haplotype_comparisons.tsv"))
  frequencies <- read_tsv(file.path(run_dir, "haplotype_group_frequencies.tsv"))
  rows <- contrast_rows(leads, comparisons, frequencies, summary)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write_rows(file.path(out, "contrast_data.tsv"), rows, CONTRAST_FIELDS)
  for (type in c("png", "pdf", "svg")) draw_haplotype_contrasts(rows, summary, file.path(out, paste0("haplotype_contrasts.", type)), type)
  ntested <- sum(vapply(rows, function(r) r$comparison_status == "tested", logical(1)))
  interpretation <- c("Disease-versus-control haplotype comparison", "",
    sprintf("Diabetes versus Control and Prediabetes versus Control: %d estimable contrasts.", ntested),
    if (!ntested) "No haplotype can be identified as different between these groups from the available validated inputs." else
      "Estimates compare participant-carrier frequencies; two chromosome copies are not independent participants.",
    "Not estimable means unavailable, not zero difference, equal frequencies or no association.",
    "Causal variants remain unresolved. Neither a cohort association nor an external GWAS lead proves causality.",
    "Positive carrier-frequency differences mean a higher frequency in the disease group than in controls.",
    "Intervals are 95% Newcombe-Wilson intervals. q values retain BH correction across all tested pairwise haplotype comparisons; the figure does not recompute a selected subset.",
    "Haplotype chromosome-copy frequencies in the table are descriptive and conditional on resolved genotypes.",
    "All outputs are intended for the private project repository; do not publish participant-level data.")
  writeLines(interpretation, file.path(out, "INTERPRETATION.txt"))
  cat("Native R contrast figure: ", length(rows), " rows; ", ntested, " estimable contrasts; ", normalizePath(out), "\n", sep = "")
  invisible(rows)
}

if (sys.nframe() == 0L) {
  argv <- commandArgs(TRUE)
  if (length(argv) != 4L || !setequal(argv[c(1, 3)], c("--run-dir", "--out"))) stop("Use: Rscript code/R/plot_haplotype_contrasts.R --run-dir COMPLETED_R_RUN --out NEW_DIRECTORY")
  args <- as.list(argv[c(2, 4)]); names(args) <- gsub("-", "_", sub("^--", "", argv[c(1, 3)]))
  do.call(export_haplotype_contrasts, args)
}
