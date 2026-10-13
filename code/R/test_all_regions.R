#!/usr/bin/env Rscript
entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(entry)), "compare_haplotypes.R"))
source(file.path(dirname(normalizePath(entry)), "plot_all_regions.R"))
source(file.path(dirname(normalizePath(entry)), "replot_all_regions.R"))
assert <- function(x, message) if (!isTRUE(x)) stop(message, call. = FALSE)
lead <- list(locus_id = "L0001", chrom = "chr1", lead_pos = "100", start = "50", end = "150", status = "lead_below_min_call_rate")
groups <- list(T2D = 4, Prediabetes = 5, Control = 6)
record_rows <- list()
for (pos in c(80, 100, 120)) for (g in names(groups)) for (allele in 0:1) {
  called <- if (pos == 120) 0 else if (g == "T2D") 2 else groups[[g]]
  alt_copies <- if (called) called else 0
  r <- setNames(rep(list("0"), length(REGIONAL_FIELDS)), REGIONAL_FIELDS)
  r$locus_id <- lead$locus_id; r$chrom <- "chr1"; r$pos <- as.character(pos); r$ref <- "A"; r$alt <- "G"; r$type <- "SNP"
  r$variant_id <- paste0("chr1:", pos, ":A:G"); r$group <- g; r$n_group <- as.character(groups[[g]])
  r$distance_to_lead_bp <- as.character(pos-as.numeric(lead$lead_pos)); r$allele <- if (allele) "G" else "A"
  r$n_called <- as.character(called); r$called_chromosomes <- as.character(2*called); r$allele_index <- as.character(allele)
  r$allele_copies <- as.character(if (allele) alt_copies else 2*called-alt_copies)
  r$is_exact_lead_match <- as.character(as.integer(pos == 100))
  record_rows[[length(record_rows)+1L]] <- r
}
records <- do.call(rbind, lapply(record_rows, as.data.frame, stringsAsFactors = FALSE))
data <- regional_plot_data(records, lead, groups)
assert(length(data$rows) == 9L, "Failed lead suppressed regional SNPs")
assert(all(vapply(data$rows, function(r) r$pos %in% c(80, 100, 120), logical(1))), "Nonlead SNPs omitted")
zero <- Filter(function(r) r$pos == 120, data$rows)
assert(all(vapply(zero, function(r) r$completeness == 0 && is.na(r$observed_alt_fraction), logical(1))), "Missing calls became ALT frequency zero or reference calls")
assert(all(vapply(Filter(function(r) r$pos != 120, data$rows), function(r) r$observed_alt_fraction == .5 && r$called_chromosomes == 2*r$n_called, logical(1))), "Wrong ALT denominator")
assert(data$interpretation[[1]]$median_completeness == .5 && data$interpretation[[1]]$sites_without_calls == 1, "Completeness summary drops missing sites")
assert(data$interpretation[[1]]$median_observed_alt_fraction == .5, "ALT median includes uncalled sites")
bad <- records; bad$allele_copies[1] <- "99"
assert(inherits(try(regional_plot_data(bad, lead, groups), silent = TRUE), "try-error"), "REF+ALT conservation not checked")
bad <- records; bad$allele[1] <- "G"
assert(inherits(try(regional_plot_data(bad, lead, groups), silent = TRUE), "try-error"), "Allele/index disagreement accepted")
bad <- records; bad$distance_to_lead_bp[1] <- "99"
assert(inherits(try(regional_plot_data(bad, lead, groups), silent = TRUE), "try-error"), "Wrong lead distance accepted")
bad <- records; bad$is_exact_lead_match[1] <- "1"
assert(inherits(try(regional_plot_data(bad, lead, groups), silent = TRUE), "try-error"), "Inconsistent lead flags accepted")
not_lead <- records; not_lead$is_exact_lead_match <- "0"
assert(is.na(regional_plot_data(not_lead, lead, groups)$interpretation[[1]]$lead_completeness), "Substitute lead silently selected")
multi <- records; multi$alt <- "G,T"
excluded <- regional_plot_data(multi, lead, groups)
assert(length(excluded$rows) == 0 && excluded$interpretation[[1]]$other_records_omitted == 3, "Multiallelic SNPs treated as biallelic")
folder <- tempfile("native_R_allregions_test_"); dir.create(folder)
run <- file.path(folder, "run"); dir.create(run)
lead2 <- lead; lead2$locus_id <- "L0002"; lead2$lead_pos <- "300"; lead2$start <- "250"; lead2$end <- "350"
write_rows(file.path(run, "lead_loci.tsv"), list(lead, lead2), names(lead))
write_rows(file.path(run, "regional_variant_frequencies.tsv.gz"), record_rows, REGIONAL_FIELDS)
summary <- list(group_counts = groups, selected_loci = 2L, distance_clusters_available = 2L, synthetic = TRUE,
                parameters = list(min_call_rate = .9), implementation = "independent_R", inputs_stable_during_read = TRUE)
jsonlite::write_json(summary, file.path(run, "run_summary.json"), auto_unbox = TRUE)
writeLines("Synthetic complete", file.path(run, "COMPLETED.txt"))
manifest <- render_all_regions(run, file.path(folder, "figures"))
assert(manifest$selected_loci == 2 && length(manifest$regions) == 2, "Empty-SNP region omitted")
assert(file.exists(file.path(folder, "figures", "L0002_chr1_300", "region.png")), "No-SNP region lacks separate picture")
redrawn <- replot_all_regions(run, file.path(folder, "figures"), file.path(folder, "redrawn"), workers = 2L)
assert(isTRUE(redrawn$saved_plot_data_unchanged) && length(redrawn$regions) == 2L, "Saved-table redraw lost a region")
original_data <- file.path(folder, "figures", "L0001_chr1_100", "plot_data.tsv.gz")
redrawn_data <- file.path(folder, "redrawn", "L0001_chr1_100", "plot_data.tsv.gz")
assert(identical(unname(tools::md5sum(original_data)), unname(tools::md5sum(redrawn_data))), "Graphics-only redraw changed numerical data")
summary$distance_clusters_available <- 3L
jsonlite::write_json(summary, file.path(run, "run_summary.json"), auto_unbox = TRUE)
assert(inherits(try(render_all_regions(run, file.path(folder, "wrong_scope")), silent = TRUE), "try-error"), "Subset accepted as all regions")
assert(inherits(try(replot_all_regions(run, file.path(folder, "figures"), file.path(folder, "wrong_redraw")), silent = TRUE), "try-error"), "Redraw accepted changed source analysis")
unlink(folder, recursive = TRUE)
cat("PASS: native R all-region scope, all SNPs, missingness, copy conservation, exact leads and no-SNP rendering.\n")
