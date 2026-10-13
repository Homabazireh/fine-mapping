#!/usr/bin/env Rscript
# Independent R rendering of every selected GWAS cluster from aggregate tables.
# Includes all available biallelic PASS SNPs, irrespective of lead eligibility.
entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
if (!exists("read_tsv", mode = "function")) source(file.path(dirname(normalizePath(entry)), "compare_haplotypes.R"))

REGION_PLOT_FIELDS <- strsplit("locus_id variant_id chrom pos ref alt distance_kb is_exact_lead_match group n_group n_called completeness alt_copies called_chromosomes observed_alt_fraction", " ")[[1]]
REGION_INTERPRETATION_FIELDS <- strsplit("locus_id chrom lead_pos start end haplotype_panel_status group n_group biallelic_snps other_records_omitted sites_with_calls sites_without_calls median_completeness median_observed_alt_fraction lead_n_called lead_completeness lead_observed_alt_fraction interpretation", " ")[[1]]

regional_plot_data <- function(records, lead, group_counts) {
  if (!all(REGIONAL_FIELDS %in% names(records))) stop("Regional table is missing required columns.")
  if (nrow(records) && any(records$locus_id != lead$locus_id)) stop("Locus contains records assigned elsewhere.")
  groups <- names(group_counts)
  if (nrow(records) && any(!records$group %in% groups)) stop("Unknown regional group.")
  is_snp <- records$type == "SNP" & nchar(records$ref) == 1L & nchar(records$alt) == 1L &
    grepl("^[ACGT]$", records$ref) & grepl("^[ACGT]$", records$alt)
  omitted <- length(unique(records$variant_id[!is_snp]))
  excluded_records <- unique(records[!is_snp, c("variant_id", "type"), drop = FALSE])
  omitted_types <- excluded_records$type
  omitted_types[omitted_types == "SNP"] <- "multiallelic_SNP"
  omitted_types <- as.list(table(omitted_types))
  sites <- records[is_snp, , drop = FALSE]
  rows <- list()
  if (nrow(sites)) {
    for (field in c("pos", "distance_to_lead_bp", "n_group", "n_called", "allele_index", "allele_copies", "called_chromosomes", "is_exact_lead_match")) {
      sites[[field]] <- suppressWarnings(as.numeric(sites[[field]]))
      if (any(!is.finite(sites[[field]])) || any(sites[[field]] != floor(sites[[field]]))) stop("Invalid integer SNP field: ", field)
    }
    if (any(sites$pos < as.numeric(lead$start) | sites$pos > as.numeric(lead$end)) || any(sites$chrom != lead$chrom)) stop("SNP outside the declared region.")
    if (any(sites$distance_to_lead_bp != sites$pos-as.numeric(lead$lead_pos))) stop("SNP distance disagrees with the selected lead.")
    variant_ids <- unique(sites$variant_id)
    by_variant <- split(seq_len(nrow(sites)), sites$variant_id)
    for (variant in variant_ids) {
      record <- sites[by_variant[[variant]], , drop = FALSE]
      if (length(unique(paste(record$chrom, record$pos, record$ref, record$alt))) != 1L) stop("Conflicting variant coordinates or alleles.")
      if (length(unique(record$is_exact_lead_match)) != 1L || !record$is_exact_lead_match[1] %in% c(0, 1)) stop("Lead flag must agree across all groups.")
      if (record$is_exact_lead_match[1] == 1 && record$pos[1] != as.numeric(lead$lead_pos)) stop("Exact lead flag at the wrong coordinate.")
      for (g in groups) {
        pair <- record[record$group == g, , drop = FALSE]
        if (nrow(pair) != 2L || !identical(sort(pair$allele_index), c(0, 1))) stop("Every SNP/group needs exactly one REF and one ALT row.")
        if (any(pair$allele != ifelse(pair$allele_index == 0, pair$ref, pair$alt))) stop("Allele index does not match the REF/ALT sequence.")
        if (length(unique(pair$n_group)) != 1L || pair$n_group[1] != group_counts[[g]] ||
            length(unique(pair$n_called)) != 1L || length(unique(pair$is_exact_lead_match)) != 1L) stop("Inconsistent SNP denominators or lead flags.")
        total <- pair$n_group[1]; called <- pair$n_called[1]
        if (called < 0 || called > total || total < 1 || any(pair$allele_copies < 0) ||
            any(pair$called_chromosomes != 2*called) || sum(pair$allele_copies) != 2*called) stop("REF and ALT copies must sum to twice the called participants.")
        alt <- pair[pair$allele_index == 1, , drop = FALSE]
        rows[[length(rows)+1L]] <- list(locus_id = lead$locus_id, variant_id = variant, chrom = alt$chrom,
          pos = alt$pos, ref = alt$ref, alt = alt$alt, distance_kb = (alt$pos-as.numeric(lead$lead_pos))/1000,
          is_exact_lead_match = alt$is_exact_lead_match, group = g, n_group = total, n_called = called,
          completeness = called/total, alt_copies = alt$allele_copies, called_chromosomes = 2*called,
          observed_alt_fraction = if (called) alt$allele_copies/(2*called) else NA_real_)
      }
    }
    order_rows <- order(vapply(rows, `[[`, numeric(1), "pos"), vapply(rows, `[[`, character(1), "variant_id"),
                        match(vapply(rows, `[[`, character(1), "group"), groups))
    rows <- rows[order_rows]
  }
  stats <- lapply(groups, function(g) {
    selected <- Filter(function(r) r$group == g, rows)
    called <- vapply(selected, `[[`, numeric(1), "n_called")
    completeness <- vapply(selected, `[[`, numeric(1), "completeness")
    alt <- vapply(selected, `[[`, numeric(1), "observed_alt_fraction")
    matched <- Filter(function(r) r$is_exact_lead_match == 1, selected)
    exact <- if (length(matched) == 1L) matched[[1]] else NULL
    nsites <- length(selected)
    interpretation <- if (!nsites) "No biallelic SNP records are available in this window; regional frequencies are unavailable." else
      paste0("Completeness is called participants divided by group size. ALT fractions condition on called chromosomes; missing calls are not reference calls. ",
             "These descriptive SNP patterns do not establish haplotype differences, cohort associations or causal variants.")
    list(locus_id = lead$locus_id, chrom = lead$chrom, lead_pos = lead$lead_pos, start = lead$start, end = lead$end,
      haplotype_panel_status = lead$status, group = g, n_group = group_counts[[g]], biallelic_snps = nsites,
      other_records_omitted = omitted, sites_with_calls = sum(called > 0), sites_without_calls = sum(called == 0),
      median_completeness = if (nsites) median(completeness) else NA_real_,
      median_observed_alt_fraction = if (any(is.finite(alt))) median(alt[is.finite(alt)]) else NA_real_,
      lead_n_called = if (is.null(exact)) NA_real_ else exact$n_called,
      lead_completeness = if (is.null(exact)) NA_real_ else exact$completeness,
      lead_observed_alt_fraction = if (is.null(exact)) NA_real_ else exact$observed_alt_fraction,
      interpretation = interpretation)
  })
  list(rows = rows, interpretation = stats, omitted_record_types = omitted_types)
}

draw_all_region <- function(data, lead, summary, path, type = "png") {
  groups <- names(summary$group_counts)
  groups <- groups[order(match(display_group(groups), c("Diabetes", "Prediabetes", "Control")), groups)]
  colors <- c(Diabetes = "#c66a32", Prediabetes = "#218d88", Control = "#3569a1")
  if (type == "png") png(path, width = 12, height = 9.5, units = "in", res = 300, type = "cairo")
  else if (type == "pdf") pdf(path, width = 12, height = 9.5, useDingbats = FALSE)
  else stop("Unsupported regional image format.")
  on.exit(dev.off())
  par(mfrow = c(length(groups), 2), mar = c(3.2, 4.8, 2.2, 1), oma = c(6, .5, 6, .5), family = "sans",
      mgp = c(2.1, .65, 0), col.axis = "#40566d", col.lab = "#233445", fg = "#233445")
  xlim <- (c(as.numeric(lead$start), as.numeric(lead$end))-as.numeric(lead$lead_pos))/1000
  if (diff(xlim) == 0) xlim <- xlim + c(-.5, .5)
  for (g in groups) {
    rows <- Filter(function(r) r$group == g, data$rows)
    x <- vapply(rows, `[[`, numeric(1), "distance_kb")
    color <- unname(colors[display_group(g)]); if (is.na(color)) color <- "#64748b"
    for (metric in c("completeness", "observed_alt_fraction")) {
      y <- vapply(rows, `[[`, numeric(1), metric)*100
      label <- if (metric == "completeness") "Genotype completeness (%)" else "ALT fraction among called chromosomes (%)"
      plot(NA_real_, xlim = xlim, ylim = c(0, 100), xlab = "Distance from selected lead (kb)", ylab = label,
           yaxs = "i", axes = FALSE, main = paste0(display_group(g), " (n = ", summary$group_counts[[g]], ")"), cex.main = 1, cex.lab = .78)
      abline(h = seq(0, 100, 25), col = "#e2e8ee", lwd = .7)
      abline(v = 0, col = "#9ba9b7", lty = 3)
      axis(1, cex.axis = .8); axis(2, at = seq(0, 100, 25), las = 1, cex.axis = .8)
      box(bty = "l", col = "#c8d1db")
      if (length(x)) {
        valid <- is.finite(y)
        points(x[valid], y[valid], pch = 16, cex = .28, col = adjustcolor(color, alpha.f = .6))
        exact <- vapply(rows, function(r) r$is_exact_lead_match == 1, logical(1)) & valid
        if (sum(exact) == 1) points(x[exact], y[exact], pch = 23, bg = color, col = "#233445", cex = 1.1)
        mtext(sprintf("%d/%d SNPs have a plotted estimate", sum(valid), length(x)), side = 3, adj = 1, line = .1, cex = .62, col = "#63758a")
      } else text(mean(xlim), 50, "No biallelic SNP records\nin this window", cex = .9, col = "#63758a")
    }
  }
  main <- paste0(lead$locus_id, "  ", lead$chrom, ":", format(as.numeric(lead$start), big.mark = ",", scientific = FALSE, trim = TRUE),
                 "-", format(as.numeric(lead$end), big.mark = ",", scientific = FALSE, trim = TRUE), " | All available regional SNPs")
  if (isTRUE(summary$synthetic)) main <- paste("SYNTHETIC FIXTURE |", main)
  mtext(main, side = 3, outer = TRUE, line = 3.6, adj = .03, font = 2, cex = 1.18, col = "#172f46")
  mtext(paste0("Lead ", lead$chrom, ":", lead$lead_pos, " | Panel status: ", gsub("_", " ", lead$status, fixed = TRUE)),
        side = 3, outer = TRUE, line = 2, adj = .03, cex = .78, col = "#526579")
  mtext("Every available PASS biallelic SNP is included; regional plots do not require an eligible lead haplotype panel.",
        side = 3, outer = TRUE, line = .5, adj = .03, cex = .73, col = "#526579")
  mtext("Dots: individual SNPs. Diamond: unique exact lead, when available. Vertical line: selected lead position.",
        side = 1, outer = TRUE, line = 1.1, adj = .03, cex = .73, col = "#526579")
  mtext("ALT fractions are conditional on called chromosomes; missing calls remain missing and can bias apparent group differences.",
        side = 1, outer = TRUE, line = 2.6, adj = .03, cex = .73, col = "#526579")
  mtext("Descriptive SNP patterns do not identify haplotype associations or causal variants. Non-SNP/multiallelic records are omitted.",
        side = 1, outer = TRUE, line = 4.1, adj = .03, cex = .73, col = "#526579")
}

split_regional_stream <- function(path, leads, folder) {
  dir.create(folder, recursive = TRUE)
  con <- gzfile(path, "rt"); on.exit(close(con))
  header <- readLines(con, n = 1L, warn = FALSE)
  if (!length(header)) stop("Empty regional table.")
  allowed <- leads$locus_id
  repeat {
    lines <- readLines(con, n = 30000L, warn = FALSE)
    if (!length(lines)) break
    text <- paste(c(header, lines), collapse = "\n")
    if (requireNamespace("data.table", quietly = TRUE)) {
      chunk <- as.data.frame(data.table::fread(text = text, sep = "\t", quote = "", colClasses = "character", na.strings = NULL, showProgress = FALSE))
    } else {
      input <- textConnection(text)
      chunk <- read.delim(input, quote = "", comment.char = "", check.names = FALSE, colClasses = "character", na.strings = character())
      close(input)
    }
    if (any(!chunk$locus_id %in% allowed)) stop("Regional records contain an unselected locus.")
    for (id in unique(chunk$locus_id)) {
      file <- file.path(folder, paste0(id, ".tsv")); exists <- file.exists(file)
      write.table(chunk[chunk$locus_id == id, , drop = FALSE], file, sep = "\t", quote = FALSE,
                  row.names = FALSE, col.names = !exists, append = exists, na = "NA")
    }
  }
}

render_all_regions <- function(run_dir, out) {
  if (!file.exists(file.path(run_dir, "COMPLETED.txt"))) stop("A completed all-cluster analysis is required.")
  if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) stop("Use a new or empty regional output directory.")
  input_paths <- file.path(run_dir, c("run_summary.json", "lead_loci.tsv", "regional_variant_frequencies.tsv.gz"))
  before <- tools::md5sum(input_paths)
  summary <- jsonlite::fromJSON(file.path(run_dir, "run_summary.json"), simplifyVector = FALSE)
  if (!isTRUE(summary$inputs_stable_during_read)) stop("The source analysis did not verify stable inputs.")
  leads <- read_tsv(file.path(run_dir, "lead_loci.tsv"))
  if (anyDuplicated(leads$locus_id) || any(!grepl("^L[0-9]+$", leads$locus_id))) stop("Invalid or duplicate locus IDs.")
  if (nrow(leads) != summary$selected_loci || is.null(summary$distance_clusters_available) || nrow(leads) != summary$distance_clusters_available)
    stop("All-region output requires every available GWAS cluster; rerun the analysis with --top-loci 0.")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  spool <- tempfile("R_regional_spool_"); on.exit(unlink(spool, recursive = TRUE), add = TRUE)
  source <- file.path(run_dir, "regional_variant_frequencies.tsv.gz")
  split_regional_stream(source, leads, spool)
  if (!file.copy(source, file.path(out, "all_regional_variant_frequencies.tsv.gz"))) stop("Could not preserve the full aggregate variant table.")
  interpretations <- list(); manifest_rows <- list()
  for (i in seq_len(nrow(leads))) {
    lead <- as.list(leads[i, ]); path <- file.path(spool, paste0(lead$locus_id, ".tsv"))
    records <- if (file.exists(path)) read_tsv(path) else as.data.frame(setNames(rep(list(character()), length(REGIONAL_FIELDS)), REGIONAL_FIELDS))
    data <- regional_plot_data(records, lead, summary$group_counts)
    name <- paste(lead$locus_id, lead$chrom, lead$lead_pos, sep = "_")
    if (!grepl("^[A-Za-z0-9_.-]+$", name)) stop("Unsafe region folder name.")
    folder <- file.path(out, name); dir.create(folder)
    write_rows(file.path(folder, "plot_data.tsv.gz"), data$rows, REGION_PLOT_FIELDS)
    for (type in c("png", "pdf")) draw_all_region(data, lead, summary, file.path(folder, paste0("region.", type)), type)
    text <- c(paste0("# ", lead$locus_id, ": ", lead$chrom, ":", lead$start, "-", lead$end), "",
      paste("Haplotype panel status:", lead$status), "", "![Regional SNP completeness and observed ALT fractions](region.png)", "",
      "All available PASS biallelic SNPs in the window are shown, including regions whose lead fails haplotype QC.", "",
      vapply(data$interpretation, function(r) {
        percent <- function(value) if (is.na(value)) "unavailable" else sprintf("%.1f%%", value*100)
        paste0(display_group(r$group), ": ", r$biallelic_snps, " SNPs; ", r$sites_with_calls,
          " have called genotypes; ", r$sites_without_calls, " have no called genotypes. Median completeness: ", percent(r$median_completeness),
          "; median observed ALT fraction among sites with calls: ", percent(r$median_observed_alt_fraction), ". ", r$interpretation)
      }, character(1)), "", paste("Omitted plotted record types:", json(data$omitted_record_types)),
      "All omitted record types remain in the complete aggregate variant table linked from the global index.", "",
      "A high observed ALT fraction can reflect missing reference calls; it is not evidence that the cohort carries only alternate alleles.",
      "These plots do not test group associations and cannot identify causal variants. Subject-level phase and independent evidence remain necessary.")
    writeLines(text, file.path(folder, "INTERPRETATION.md"))
    interpretations <- c(interpretations, data$interpretation)
    manifest_rows[[i]] <- list(locus_id = lead$locus_id, chrom = lead$chrom, lead_pos = lead$lead_pos,
      directory = name, biallelic_snps = data$interpretation[[1]]$biallelic_snps, other_records_omitted = data$interpretation[[1]]$other_records_omitted,
      omitted_record_types = data$omitted_record_types,
      files = list(png = paste0(name, "/region.png"), pdf = paste0(name, "/region.pdf"), data = paste0(name, "/plot_data.tsv.gz"), interpretation = paste0(name, "/INTERPRETATION.md")))
    if (i %% 25L == 0L || i == nrow(leads)) cat("Independent R regional figures: ", i, "/", nrow(leads), "\n", sep = "")
  }
  write_rows(file.path(out, "interpretations.tsv"), interpretations, REGION_INTERPRETATION_FIELDS)
  index <- c("# All regional SNP figures", "", "Private project outputs. Each region has a separate PNG, PDF, full plotting table and interpretation.", "",
    "All GWAS clusters are represented, including those with unavailable haplotypes or no available biallelic SNPs.", "",
    "ALT fractions use only called chromosomes; no missing genotype is converted to a reference call. These pictures do not establish haplotype or causal associations.", "",
    "[Full aggregate variant table, including omitted indels and multiallelic sites](all_regional_variant_frequencies.tsv.gz)", "",
    "| Region | Window lead | SNPs | Interpretation | PNG | PDF |", "| --- | --- | ---: | --- | --- | --- |",
    vapply(manifest_rows, function(r) paste0("| ", r$locus_id, " | ", r$chrom, ":", r$lead_pos, " | ", r$biallelic_snps,
      " | [Read](", r$files$interpretation, ") | [PNG](", r$files$png, ") | [PDF](", r$files$pdf, ") |"), character(1)))
  writeLines(index, file.path(out, "index.md"))
  html <- c("<!doctype html><html lang='en'><meta charset='utf-8'><title>Private regional SNP figures</title>",
    "<style>body{font:16px sans-serif;max-width:1100px;margin:3rem auto;color:#233445}table{border-collapse:collapse;width:100%}td,th{padding:.5rem;border-bottom:1px solid #ddd;text-align:left}</style>",
    "<h1>All regional SNP figures</h1><p>Private outputs. Each region contains all available biallelic SNPs, irrespective of lead eligibility. ALT fractions condition on called chromosomes; these descriptive pictures do not establish haplotype associations or causal variants.</p>",
    "<p><a href='all_regional_variant_frequencies.tsv.gz'>Full aggregate variant table, including omitted indels and multiallelic records</a></p>",
    "<table><tr><th>Region</th><th>Lead</th><th>SNPs</th><th>Files</th></tr>",
    vapply(manifest_rows, function(r) paste0("<tr><td>", r$locus_id, "</td><td>", r$chrom, ":", r$lead_pos, "</td><td>", r$biallelic_snps,
      "</td><td><a href='", r$files$png, "'>PNG</a> &middot; <a href='", r$files$pdf, "'>PDF</a> &middot; <a href='", r$files$interpretation,
      "'>Interpretation</a> &middot; <a href='", r$files$data, "'>Data</a></td></tr>"), character(1)), "</table></html>")
  writeLines(html, file.path(out, "index.html"))
  if (!identical(before, tools::md5sum(input_paths))) stop("An aggregate input changed during plotting; discard these outputs.")
  manifest <- list(implementation = "independent_R", created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    synthetic = isTRUE(summary$synthetic), selected_loci = nrow(leads), distance_clusters_available = summary$distance_clusters_available,
    all_clusters_included = TRUE, inputs_stable_during_plotting = TRUE, png_dpi = 300, estimates = "SNP completeness and ALT fraction conditional on called chromosomes; no inference",
    full_regional_variant_table = "all_regional_variant_frequencies.tsv.gz", input_md5 = as.list(setNames(unname(before), basename(input_paths))),
    input_source = if (identical(summary$implementation, "independent_R")) "independent_R_aggregate_run" else "Python_aggregate_run_independent_R_validation_and_plotting",
    versions = list(R = R.version.string, jsonlite = as.character(packageVersion("jsonlite"))), regions = manifest_rows)
  jsonlite::write_json(manifest, file.path(out, "manifest.json"), pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null")
  invisible(manifest)
}

if (sys.nframe() == 0L) {
  argv <- commandArgs(TRUE)
  if (length(argv) != 4L || !setequal(argv[c(1, 3)], c("--run-dir", "--out"))) stop("Use: Rscript code/R/plot_all_regions.R --run-dir COMPLETED_ALL_CLUSTER_RUN --out NEW_DIRECTORY")
  args <- as.list(argv[c(2, 4)]); names(args) <- gsub("-", "_", sub("^--", "", argv[c(1, 3)]))
  do.call(render_all_regions, args)
}
