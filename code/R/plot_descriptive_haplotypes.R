#!/usr/bin/env Rscript
# Observed chromosome-copy frequencies only. No association or causal inference.
entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
if (!exists("read_tsv", mode = "function")) source(file.path(dirname(normalizePath(entry)), "compare_haplotypes.R"))

DESCRIPTIVE_PLOT_FIELDS <- strsplit("locus_id chrom lead_pos haplotype_id group n_group n_resolved haplotype_copies resolved_chromosomes haplotype_frequency", " ")[[1]]
DESCRIPTIVE_DIFFERENCE_FIELDS <- strsplit("locus_id chrom lead_pos haplotype_id contrast case_group control_group haplotype_copies_case haplotype_copies_control resolved_chromosomes_case resolved_chromosomes_control n_resolved_case n_resolved_control n_group_case n_group_control haplotype_frequency_case haplotype_frequency_control descriptive_copy_frequency_difference interpretation_status", " ")[[1]]
DESCRIPTIVE_REGION_FIELDS <- strsplit("locus_id contrast haplotype_id_largest_absolute_descriptive_difference descriptive_copy_frequency_difference n_resolved_case n_resolved_control interpretation", " ")[[1]]

descriptive_haplotype_data <- function(leads, frequencies, summary) {
  groups <- names(summary$group_counts)
  counts <- unlist(summary$group_counts)
  if (length(groups) != 3L || anyDuplicated(groups) || any(!is.finite(counts) | counts <= 0 | counts != floor(counts)))
    stop("Descriptive figure requires exactly three distinct groups with positive integer counts.")
  if (anyDuplicated(leads$locus_id)) stop("Duplicate lead IDs.")
  lookup <- function(label) {
    matches <- groups[display_group(groups) == label]
    if (length(matches) != 1L) stop("Expected one explicit source group for ", label, "; groups are never pooled.")
    matches
  }
  groups <- setNames(vapply(c("Diabetes", "Prediabetes", "Control"), lookup, character(1)), c("Diabetes", "Prediabetes", "Control"))
  required <- c("locus_id", "haplotype_id", "group", "n_group", "n_resolved", "haplotype_copies", "resolved_chromosomes", "haplotype_frequency")
  if (!all(required %in% names(frequencies))) stop("Incomplete haplotype frequency schema.")
  if (!nrow(frequencies)) stop("No resolved haplotype frequencies are available.")
  if (anyDuplicated(paste(frequencies$locus_id, frequencies$haplotype_id, frequencies$group))) stop("Duplicate locus/haplotype/group frequency.")
  if (any(!frequencies$locus_id %in% leads$locus_id) || any(!frequencies$group %in% groups)) stop("Frequency references an unknown region or group.")
  loci <- leads$locus_id[leads$locus_id %in% frequencies$locus_id]
  if (any(leads$status[match(loci, leads$locus_id)] != "analysed")) stop("Haplotype frequencies supplied for an ineligible panel.")
  data <- frequencies
  for (field in c("n_group", "n_resolved", "haplotype_copies", "resolved_chromosomes")) {
    data[[field]] <- suppressWarnings(as.numeric(data[[field]]))
    if (any(!is.finite(data[[field]])) || any(data[[field]] < 0 | data[[field]] != floor(data[[field]]))) stop("Invalid count: ", field)
  }
  if (any(data$n_group < 1) || any(data$n_resolved > data$n_group) || any(data$resolved_chromosomes != 2*data$n_resolved) || any(data$haplotype_copies > data$resolved_chromosomes)) stop("Invalid resolved participant/chromosome denominators.")
  provided <- suppressWarnings(as.numeric(data$haplotype_frequency))
  expected <- ifelse(data$resolved_chromosomes > 0, data$haplotype_copies/data$resolved_chromosomes, NA_real_)
  if (!identical(is.na(provided), is.na(expected)) || any(abs(provided[is.finite(expected)]-expected[is.finite(expected)]) > 1e-7)) stop("Stored haplotype frequency disagrees with copy counts.")
  data$haplotype_frequency <- expected
  plot <- differences <- interpretation <- list()
  for (locus in loci) {
    lead <- as.list(leads[match(locus, leads$locus_id), ])
    rows <- data[data$locus_id == locus, , drop = FALSE]
    haps <- sort(unique(rows$haplotype_id))
    for (g in unname(groups)) {
      subset <- rows[rows$group == g, , drop = FALSE]
      if (nrow(subset) != length(haps) || !setequal(subset$haplotype_id, haps)) stop("Every local haplotype needs a row for every group.")
      if (length(unique(subset$n_group)) != 1L || unique(subset$n_group) != summary$group_counts[[g]] ||
          length(unique(subset$n_resolved)) != 1L || length(unique(subset$resolved_chromosomes)) != 1L) stop("Inconsistent group/phase denominators across local haplotypes.")
      if (sum(subset$haplotype_copies) != subset$resolved_chromosomes[1]) stop("Haplotype copies do not sum to resolved chromosomes.")
    }
    for (hap in haps) {
      select <- function(g) as.list(rows[rows$haplotype_id == hap & rows$group == g, , drop = FALSE][1, ])
      for (g in unname(groups)) {
        value <- select(g)
        plot[[length(plot)+1L]] <- list(locus_id = locus, chrom = lead$chrom, lead_pos = lead$lead_pos, haplotype_id = hap,
          group = g, n_group = value$n_group, n_resolved = value$n_resolved, haplotype_copies = value$haplotype_copies,
          resolved_chromosomes = value$resolved_chromosomes, haplotype_frequency = value$haplotype_frequency)
      }
      control <- select(groups[["Control"]])
      for (label in c("Diabetes", "Prediabetes")) {
        case <- select(groups[[label]])
        differences[[length(differences)+1L]] <- list(locus_id = locus, chrom = lead$chrom, lead_pos = lead$lead_pos, haplotype_id = hap,
          contrast = paste(label, "vs Control"), case_group = groups[[label]], control_group = groups[["Control"]],
          haplotype_copies_case = case$haplotype_copies, haplotype_copies_control = control$haplotype_copies,
          resolved_chromosomes_case = case$resolved_chromosomes, resolved_chromosomes_control = control$resolved_chromosomes,
          n_resolved_case = case$n_resolved, n_resolved_control = control$n_resolved, n_group_case = case$n_group, n_group_control = control$n_group,
          haplotype_frequency_case = case$haplotype_frequency, haplotype_frequency_control = control$haplotype_frequency,
          descriptive_copy_frequency_difference = case$haplotype_frequency-control$haplotype_frequency,
          interpretation_status = "descriptive_only_not_tested")
      }
    }
    for (label in c("Diabetes", "Prediabetes")) {
      selected <- Filter(function(r) r$locus_id == locus && r$contrast == paste(label, "vs Control"), differences)
      finite <- Filter(function(r) is.finite(r$descriptive_copy_frequency_difference), selected)
      best <- NULL
      if (length(finite)) {
        size <- abs(vapply(finite, `[[`, numeric(1), "descriptive_copy_frequency_difference"))
        tied <- finite[abs(size-max(size)) <= 1e-12]
        best <- tied[[order(vapply(tied, `[[`, character(1), "haplotype_id"))[1]]]
      }
      message <- if (is.null(best)) "Descriptive contrast unavailable because one or both groups have no resolved chromosome copies." else
        sprintf("Largest absolute observed chromosome-copy frequency contrast: %s, %s vs Control, %+.2f percentage points among %d/%d resolved participants. Descriptive only; selected observed calls, no association test or causal conclusion.",
          best$haplotype_id, label, 100*best$descriptive_copy_frequency_difference, best$n_resolved_case, best$n_resolved_control)
      interpretation[[length(interpretation)+1L]] <- list(locus_id = locus, contrast = paste(label, "vs Control"),
        haplotype_id_largest_absolute_descriptive_difference = if (is.null(best)) NA_character_ else best$haplotype_id,
        descriptive_copy_frequency_difference = if (is.null(best)) NA_real_ else best$descriptive_copy_frequency_difference,
        n_resolved_case = selected[[1]]$n_resolved_case, n_resolved_control = selected[[1]]$n_resolved_control, interpretation = message)
    }
  }
  list(plot = plot, differences = differences, interpretation = interpretation, loci = loci, groups = groups)
}

draw_descriptive_haplotypes <- function(data, summary, path, type = "png", columns = 3L) {
  loci <- data$loci
  if (columns < 1 || columns != floor(columns)) stop("Facet columns must be a positive integer.")
  nrows <- max(1L, ceiling(length(loci)/columns))
  by_row <- split(seq_along(loci), ceiling(seq_along(loci)/columns))
  max_haps <- vapply(by_row, function(ix) max(vapply(loci[ix], function(locus) length(unique(vapply(Filter(function(r) r$locus_id == locus, data$plot), `[[`, character(1), "haplotype_id"))), integer(1))), integer(1))
  if (!length(max_haps)) max_haps <- 1L
  row_heights <- pmax(2.9, 1.35 + .19*max_haps)
  width <- max(10, 5.1*columns); height <- 2.35+sum(row_heights)
  if (type == "png") png(path, width = width, height = height, units = "in", res = 220, type = "cairo")
  else if (type == "pdf") pdf(path, width = width, height = height, useDingbats = FALSE)
  else if (type == "svg") svg(path, width = width, height = height)
  else stop("Unsupported descriptive figure format.")
  on.exit(dev.off())
  par(mar = rep(0, 4), family = "sans", bg = "white")
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  text(.025, 1-.22/height, "Observed haplotype frequencies and descriptive differences", adj = 0, font = 2, cex = 1.6, col = "#172f46")
  text(.025, 1-.52/height, sprintf("%d resolved regions | %d local haplotype patterns | Conditional chromosome-copy frequencies; no association tests",
       length(loci), length(data$plot)/3), adj = 0, cex = .96, col = "#526579")
  text(.025, 1-.8/height, "Gray cells: unavailable. Differences are Diabetes minus Control and Prediabetes minus Control; no CI, p-values or causal inference.",
       adj = 0, cex = .86, col = "#526579")
  sequential <- colorRampPalette(c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"))(101)
  diverging <- colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(201)
  top <- 1-1.25/height; margins <- .025; gap <- .025; facet_width <- (1-2*margins-gap*(columns-1L))/columns
  cumulative <- c(0, cumsum(row_heights))/height
  for (i in seq_along(loci)) {
    locus <- loci[i]; grid_row <- (i-1L) %/% columns+1L; grid_col <- (i-1L) %% columns+1L
    left <- margins+(grid_col-1L)*(facet_width+gap); right <- left+facet_width
    facet_top <- top-cumulative[grid_row]; facet_bottom <- top-cumulative[grid_row+1L]
    rows <- Filter(function(r) r$locus_id == locus, data$plot)
    haps <- sort(unique(vapply(rows, `[[`, character(1), "haplotype_id")))
    first <- rows[[1]]
    text(left, facet_top-.1/height, paste0(locus, "  ", first$chrom, ":", first$lead_pos), adj = 0, font = 2, cex = .94, col = "#172f46")
    denoms <- vapply(unname(data$groups), function(g) {
      r <- Filter(function(r) r$group == g, rows)[[1]]; paste0(r$n_resolved, "/", r$n_group)
    }, character(1))
    text(left, facet_top-.35/height, paste("Resolved people D/P/C:", paste(denoms, collapse = " | ")), adj = 0, cex = .7, col = "#526579")
    matrix_left <- left+.05; matrix_right <- right-.004; cell_width <- (matrix_right-matrix_left)/5
    x <- matrix_left+(seq_len(5)-.5)*cell_width
    text(x, facet_top-.64/height, c("Diabetes", "Prediabetes", "Control", "D - C", "P - C"), cex = .66, col = "#40566d")
    matrix_top <- facet_top-.84/height; matrix_bottom <- facet_bottom+.22/height
    step <- (matrix_top-matrix_bottom)/length(haps)
    for (j in seq_along(haps)) {
      hap <- haps[j]; y <- matrix_top-(j-.5)*step
      freq <- vapply(unname(data$groups), function(g) Filter(function(r) r$haplotype_id == hap && r$group == g, rows)[[1]]$haplotype_frequency, numeric(1))
      values <- c(freq, freq[1]-freq[3], freq[2]-freq[3])
      text(matrix_left-.004, y, sub(paste0(locus, "_"), "", hap, fixed = TRUE), adj = 1, cex = .63, col = "#233445")
      for (k in seq_along(values)) {
        value <- values[k]
        color <- if (is.na(value)) "#d9dde3" else if (k <= 3) sequential[round(value*100)+1L] else diverging[round((value+1)*100)+1L]
        rect(x[k]-.48*cell_width, y-.47*step, x[k]+.48*cell_width, y+.47*step, col = color, border = NA)
        label <- if (is.na(value)) "NA" else if (k <= 3) sprintf("%.2f", value) else sprintf("%+.2f", value)
        font_color <- if (!is.na(value) && ((k <= 3 && value > .65) || (k > 3 && abs(value) > .65))) "white" else "#233445"
        text(x[k], y, label, cex = .62, col = font_color)
      }
    }
  }
  if (!length(loci)) text(.5, .55, "No resolved haplotype frequencies are available.", cex = 1.4, col = "#526579")
  # Two independent color scales, clearly separated from inference.
  legend_y <- .64/height
  for (k in seq_len(101)) rect(.05+(k-1)*.002, legend_y, .05+k*.002, legend_y+.13/height, col = sequential[k], border = NA)
  text(c(.05, .15, .252), legend_y-.11/height, c("0", "Copy frequency", "1"), cex = .76, col = "#526579")
  for (k in seq_len(201)) rect(.43+(k-1)*.0012, legend_y, .43+k*.0012, legend_y+.13/height, col = diverging[k], border = NA)
  text(c(.43, .551, .671), legend_y-.11/height, c("-1", "Descriptive difference", "+1"), cex = .76, col = "#526579")
  text(.985, .25/height, paste(if (isTRUE(summary$synthetic)) "SYNTHETIC" else "Real cohort", "| Independent R | Denominator: twice resolved people | Missing/phase exclusions can bias frequencies | Causal variants unresolved"),
       adj = 1, cex = .7, col = "#63758a")
}

export_descriptive_haplotypes <- function(run_dir, out, columns = 3L) {
  if (!file.exists(file.path(run_dir, "COMPLETED.txt"))) stop("A completed analysis is required.")
  if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) stop("Use a new or empty descriptive output directory.")
  paths <- file.path(run_dir, c("run_summary.json", "lead_loci.tsv", "haplotype_group_frequencies.tsv"))
  before <- tools::md5sum(paths)
  summary <- jsonlite::fromJSON(paths[1], simplifyVector = FALSE)
  if (!isTRUE(summary$inputs_stable_during_read)) stop("Source analysis did not verify stable inputs.")
  leads <- read_tsv(paths[2])
  if (is.null(summary$distance_clusters_available) || nrow(leads) != summary$distance_clusters_available || nrow(leads) != summary$selected_loci)
    stop("The descriptive all-region figure requires every available GWAS cluster.")
  if ((summary$haplotype_carrier_tests_performed %||% 0) > 0 || (summary$haplotype_overall_tests_performed %||% 0) > 0)
    stop("This descriptive-only figure is for runs without association tests; use the tested-contrast workflow for inferential results.")
  data <- descriptive_haplotype_data(leads, read_tsv(paths[3]), summary)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write_rows(file.path(out, "plot_data.tsv"), data$plot, DESCRIPTIVE_PLOT_FIELDS)
  write_rows(file.path(out, "descriptive_frequency_differences.tsv"), data$differences, DESCRIPTIVE_DIFFERENCE_FIELDS)
  write_rows(file.path(out, "region_interpretations.tsv"), data$interpretation, DESCRIPTIVE_REGION_FIELDS)
  for (type in c("png", "pdf", "svg")) draw_descriptive_haplotypes(data, summary, file.path(out, paste0("descriptive_haplotypes.", type)), type, columns)
  writeLines(c("Observed haplotype frequency differences: descriptive only", "",
    sprintf("%d regions, %d local haplotype patterns, %d disease-control descriptive differences.", length(data$loci), length(data$plot)/3, length(data$differences)),
    "Frequencies use observed haplotype copies divided by twice the number of resolved participants in that region/group.",
    "The two difference columns subtract the Control copy frequency from Diabetes and Prediabetes, respectively.",
    "These are descriptive differences among selected, callable, phase-resolved participants. They are not tested associations or population estimates.",
    "No confidence interval, p-value, q-value or significance claim is attached to these descriptive differences.",
    "Missingness, unconfirmed reference-call provenance and phase resolution can bias both frequencies and differences.",
    "The largest absolute difference per region is reported for description only; it does not identify a significant or causal haplotype.",
    "Causal variants remain unresolved; external GWAS association does not establish causality.",
    "All outputs belong in the private project repository."), file.path(out, "INTERPRETATION.txt"))
  if (!identical(before, tools::md5sum(paths))) stop("An input changed during descriptive plotting; discard these outputs.")
  manifest <- list(implementation = "independent_R", synthetic = isTRUE(summary$synthetic),
    created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), regions = length(data$loci), haplotype_patterns = length(data$plot)/3,
    plotted_group_rows = length(data$plot), descriptive_contrasts = length(data$differences), inference_performed = FALSE,
    interpretation = "Descriptive chromosome-copy frequencies only; no association or causal inference", png_dpi = 220,
    source_implementation = summary$implementation %||% "Python", inputs_stable_during_plotting = TRUE,
    input_md5 = as.list(setNames(unname(before), basename(paths))))
  jsonlite::write_json(manifest, file.path(out, "manifest.json"), pretty = TRUE, auto_unbox = TRUE)
  cat("Native R descriptive figure: ", length(data$loci), " regions, ", length(data$plot)/3, " patterns; no inference.\n", sep = "")
  invisible(data)
}

if (sys.nframe() == 0L) {
  argv <- commandArgs(TRUE)
  if (!length(argv) || length(argv) %% 2L) stop("Use --run-dir COMPLETED_RUN --out NEW_DIRECTORY [--columns 3]")
  flags <- argv[seq(1, length(argv), 2)]; args <- as.list(argv[seq(2, length(argv), 2)])
  if (!all(flags %in% c("--run-dir", "--out", "--columns")) || anyDuplicated(flags)) stop("Unknown or repeated descriptive option.")
  names(args) <- gsub("-", "_", sub("^--", "", flags))
  if (!is.null(args$columns)) args$columns <- as.integer(args$columns)
  do.call(export_descriptive_haplotypes, args)
}
