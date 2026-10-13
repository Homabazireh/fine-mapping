#!/usr/bin/env Rscript
entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(entry)), "compare_haplotypes.R"))
source(file.path(dirname(normalizePath(entry)), "plot_haplotype_contrasts.R"))
assert <- function(x, message) if (!isTRUE(x)) stop(message, call. = FALSE)
near <- function(a, b) assert(isTRUE(all.equal(a, b, tolerance = 1e-12)), "Numeric contrast mismatch")
fails <- function(expression, message) assert(inherits(try(force(expression), silent = TRUE), "try-error"), message)
empty <- function(columns) as.data.frame(setNames(rep(list(character()), length(columns)), columns), stringsAsFactors = FALSE)
records <- function(rows, fields) do.call(rbind, lapply(rows, function(row) {
  output <- setNames(rep(list("NA"), length(fields)), fields); output[names(row)] <- row
  as.data.frame(output, stringsAsFactors = FALSE)
}))
summary_for <- function(leads, reference = TRUE) list(group_counts = list(T2D = 20, Prediabetes = 24, Control = 30),
  selected_loci = nrow(leads), locus_status_counts = as.list(table(leads$status)), inputs_stable_during_read = TRUE,
  synthetic = TRUE, parameters = list(min_call_rate = .9, test_associations = TRUE, reference_genotypes_confirmed = reference,
  min_test_samples = 10, min_resolved_fraction = .8, min_carriers = 3))
leads <- data.frame(locus_id = c("L0001", "L0002"), chrom = c("chr1", "chr2"), lead_pos = c("100", "200"),
                    status = c("lead_below_min_call_rate", "palindromic_lead_orientation_unconfirmed"))
comparisons <- empty(COMP_FIELDS); frequencies <- empty(FREQ_FIELDS); summary <- summary_for(leads)
rows <- contrast_rows(leads, comparisons, frequencies, summary)
assert(length(rows) == 4L, "Skipped loci must appear in both contrasts")
assert(all(vapply(rows, function(r) r$comparison_status == "not_estimable" && r$inference_status == "not_tested" &&
  is.na(r$descriptive_carrier_frequency_difference) && is.na(r$carrier_frequency_difference) && is.na(r$CI95_low), logical(1))), "Unavailable contrasts became zero effects")
assert(identical(vapply(rows, `[[`, character(1), "reason"), rep(leads$status, each = 2)), "Exclusion reasons changed")

lead <- leads[1, , drop = FALSE]; lead$status <- "analysed"; summary <- summary_for(lead)
comparisons <- records(list(
  list(locus_id = "L0001", haplotype_id = "L0001_H001", group_a = "Control", group_b = "T2D",
    n_resolved_a = "30", n_resolved_b = "20", carriers_a = "6", carriers_b = "8", carrier_frequency_a = "0.2", carrier_frequency_b = "0.4",
    test_status = "tested_unadjusted", carrier_frequency_difference_a_minus_b = "-0.2", carrier_difference_CI95_low = "-0.38",
    carrier_difference_CI95_high = "-0.03", p_fisher_two_sided = "0.02", q_BH_all_pairwise_tests = "0.08"),
  list(locus_id = "L0001", haplotype_id = "L0001_H001", group_a = "Prediabetes", group_b = "Control",
    n_resolved_a = "24", n_resolved_b = "30", carriers_a = "6", carriers_b = "6", carrier_frequency_a = "0.25", carrier_frequency_b = "0.2",
    test_status = "not_tested_reference_provenance")), COMP_FIELDS)
frequencies <- data.frame(locus_id = "L0001", haplotype_id = "L0001_H001", group = c("T2D", "Prediabetes", "Control"), haplotype_frequency = c("0.3", "0.2", "0.1"))
rows <- contrast_rows(lead, comparisons, frequencies, summary); disease <- rows[[1]]; prediabetes <- rows[[2]]
near(disease$carrier_frequency_difference, .2); near(disease$descriptive_carrier_frequency_difference, .2)
near(disease$CI95_low, .03); near(disease$CI95_high, .38); near(disease$q_BH_all_pairwise_tests, .08)
near(disease$n_resolved_case, 20); near(disease$carriers_case, 8); near(disease$haplotype_frequency_case_descriptive, .3)
assert(prediabetes$comparison_status == "descriptive_only" && is.na(prediabetes$carrier_frequency_difference) && is.na(prediabetes$CI95_low) && is.na(prediabetes$p_fisher_two_sided), "Descriptive observation became an inference")
near(prediabetes$descriptive_carrier_frequency_difference, .05)
assert(prediabetes$inference_status == "not_tested_reference_provenance", "Provenance reason missing")
forward <- comparisons
for (prefix in c("group", "n_resolved", "carriers", "carrier_frequency")) {
  a <- forward[[paste0(prefix, "_a")]][1]; forward[[paste0(prefix, "_a")]][1] <- forward[[paste0(prefix, "_b")]][1]; forward[[paste0(prefix, "_b")]][1] <- a
}
forward[1, c("carrier_frequency_difference_a_minus_b", "carrier_difference_CI95_low", "carrier_difference_CI95_high")] <- c("0.2", "0.03", "0.38")
near(contrast_rows(lead, forward, frequencies, summary)[[1]]$carrier_frequency_difference, disease$carrier_frequency_difference)

# Source inference is only admissible when reference, phase and carrier gates were met.
fails(contrast_rows(lead, comparisons, frequencies, summary_for(lead, FALSE)), "Provenance gate bypassed")
bad_summary <- summary; bad_summary$parameters$min_resolved_fraction <- 1.01
fails(contrast_rows(lead, comparisons, frequencies, bad_summary), "Resolution requirement bypassed")
bad_summary <- summary; bad_summary$parameters$min_carriers <- 15
fails(contrast_rows(lead, comparisons, frequencies, bad_summary), "Carrier-count requirement bypassed")
bad <- comparisons; bad$carrier_difference_CI95_high[1] <- "NA"
fails(contrast_rows(lead, bad, frequencies, summary), "Tested contrast with missing uncertainty accepted")
bad <- comparisons; bad$carrier_frequency_difference_a_minus_b[1] <- "-0.1"
fails(contrast_rows(lead, bad, frequencies, summary), "Effect disagrees with counts")
bad <- comparisons; bad$q_BH_all_pairwise_tests[1] <- "0.01"
fails(contrast_rows(lead, bad, frequencies, summary), "q smaller than p accepted")
for (field in CONTRAST_INFERENCE_FIELDS) {
  bad <- comparisons; bad[[field]][2] <- "0.01"
  fails(contrast_rows(lead, bad, frequencies, summary), paste("Untested inference value accepted:", field))
}

# Descriptive counts are validated and retain observed zero separately from missing.
descriptive <- comparisons[2, , drop = FALSE]
bad <- descriptive; bad$carriers_a <- "25"
fails(contrast_rows(lead, bad, frequencies, summary), "More carriers than participants accepted")
bad <- descriptive; bad$n_resolved_a <- "24.5"
fails(contrast_rows(lead, bad, frequencies, summary), "Fractional participant count accepted")
bad <- descriptive; bad$n_resolved_a <- "25"
fails(contrast_rows(lead, bad, frequencies, summary), "Resolved count exceeding group total accepted")
bad <- descriptive; bad$carrier_frequency_a <- "0.5"
fails(contrast_rows(lead, bad, frequencies, summary), "Inconsistent descriptive frequency accepted")
zero <- descriptive; zero$carriers_a <- "0"; zero$carriers_b <- "0"; zero$carrier_frequency_a <- "0"; zero$carrier_frequency_b <- "0"
zero_rows <- contrast_rows(lead, zero, frequencies, summary)
near(zero_rows[[2]]$descriptive_carrier_frequency_difference, 0)
assert(zero_rows[[2]]$comparison_status == "descriptive_only", "Observed zero was marked unavailable")
missing <- zero; missing$n_resolved_a <- "0"; missing$carrier_frequency_a <- "NA"
missing_rows <- contrast_rows(lead, missing, frequencies, summary)
assert(missing_rows[[2]]$comparison_status == "not_estimable" && is.na(missing_rows[[2]]$descriptive_carrier_frequency_difference), "Missing denominator became zero difference")
recomputed <- descriptive; recomputed$carrier_frequency_a <- "NA"
fails(contrast_rows(lead, recomputed, frequencies, summary), "Missing source frequency accepted")
partial <- descriptive; partial$n_resolved_a <- "NA"
fails(contrast_rows(lead, partial, frequencies, summary), "Partial participant counts accepted")
all_missing <- descriptive; all_missing[1, c("n_resolved_a", "carriers_a", "carrier_frequency_a")] <- "NA"
assert(is.na(contrast_rows(lead, all_missing, frequencies, summary)[[2]]$descriptive_carrier_frequency_difference), "Wholly missing counts became an observation")
missing_pair <- contrast_rows(lead, comparisons[0, ], frequencies, summary)
assert(all(vapply(missing_pair, function(r) r$reason == "pairwise_result_unavailable" && is.na(r$descriptive_carrier_frequency_difference), logical(1))), "Missing pair became an effect")
no_haps <- contrast_rows(lead, comparisons[0, ], frequencies[0, ], summary)
assert(all(vapply(no_haps, function(r) r$reason == "no_resolved_haplotypes", logical(1))), "Empty resolved panel reason lost")

# Interpretation distinguishes missing actual-haplotype contrasts from region
# placeholders, while retaining observed zero, unknown data and subset scope.
mixed_leads <- rbind(lead, leads[2, ], transform(lead, locus_id = "L0003", chrom = "chr3", lead_pos = "300"))
mixed_summary <- summary_for(mixed_leads, FALSE)
mixed_summary$scan_counts <- list(high_quality_homozygous_reference_genotypes = 0L)
mixed <- contrast_rows(mixed_leads, zero, frequencies, mixed_summary)
before <- serialize(mixed, NULL)
notes <- contrast_interpretation(mixed, mixed_summary)
contains <- function(pattern, text = notes) any(grepl(pattern, text, fixed = TRUE))
assert(contains("0 tested haplotype comparisons; 1 descriptive-only differences; 1 identified-haplotype comparisons"), "Haplotype missingness conflated with placeholders")
assert(contains("2 excluded-region placeholder rows cover 1 excluded regions"), "Excluded-region placeholder count wrong")
assert(contains("2 further placeholder rows cover 1 analysed regions with no resolved haplotypes"), "Unresolved analysed regions conflated with exclusions")
assert(contains("1 local haplotypes in 1 regions"), "Observed haplotype count wrong")
assert(contains("Report scope: 3 of 3 selected regions; 2 analysed regions in this report") && contains("complete source has 2 analysed (panel-eligible) regions"), "Report or source scope wrong")
assert(contains("1: pairwise_result_unavailable") && contains("1: palindromic_lead_orientation_unconfirmed") && !contains("2: palindromic_lead_orientation_unconfirmed"), "Reason count units are incorrect")
assert(contains("source analysis completed") && contains("not a plotting failure"), "Analysis status unclear")
assert(contains("reference_genotypes_confirmed=false") && contains("Source scan recorded 0 high-quality homozygous-reference genotypes"), "Source provenance evidence missing")
assert(contains("first failing gate; other limitations may also apply") && contains("participant-resolution and carrier-count requirements must also pass"), "First provenance failure hides other gates")
assert(contains("does not establish whether reference calls were absent before QC"), "Post-QC counts overinterpreted")
assert(identical(serialize(mixed, NULL), before), "Interpretation mutated input rows")
assert(is.na(mixed[[1]]$descriptive_carrier_frequency_difference) && mixed[[2]]$descriptive_carrier_frequency_difference == 0 && all(vapply(mixed[3:6], function(r) is.na(r$haplotype_id), logical(1))), "Missing, zero or unidentified haplotypes changed")
regional <- Filter(function(r) r$locus_id == "L0001", mixed)
unknown_summary <- mixed_summary; unknown_summary$parameters$reference_genotypes_confirmed <- NULL
for (count in list(NULL, NA, "NA", FALSE)) {
  unknown_summary$scan_counts <- list(high_quality_homozygous_reference_genotypes = count)
  regional_notes <- contrast_interpretation(regional, unknown_summary)
  assert(contains("Report scope: 1 of 3 selected regions; 1 analysed regions in this report", regional_notes) && contains("complete source has 2 analysed (panel-eligible) regions", regional_notes), "Regional scope replaced source scope")
  assert(contains("0 excluded-region placeholder rows cover 0 excluded regions in this report", regional_notes), "Excluded source regions incorrectly assigned to subset")
  assert(!contains("Source scan recorded", regional_notes) && !contains("reference_genotypes_confirmed=false", regional_notes), "Unknown source provenance became false or zero")
}

# Duplicate, mismatched and ineligible source records are rejected before drawing.
fails(contrast_rows(lead, rbind(comparisons, forward[1, ]), frequencies, summary), "Reversed duplicate accepted")
fails(contrast_rows(lead, comparisons, rbind(frequencies, frequencies[1, ]), summary), "Duplicate frequency accepted")
bad <- comparisons; bad$locus_id[2] <- "L9999"
fails(contrast_rows(lead, bad, frequencies, summary), "Unknown source locus accepted")
bad <- frequencies; bad$group[1] <- "Other"
fails(contrast_rows(lead, comparisons, bad, summary), "Unknown frequency group accepted")
bad_summary <- summary; bad_summary$group_counts$Diabetes <- 10
fails(contrast_rows(lead, comparisons, frequencies, bad_summary), "Ambiguous group alias accepted")
fails(contrast_rows(leads, comparisons, frequencies, summary_for(leads)), "Results attached to excluded locus accepted")
bad_summary <- summary; bad_summary$locus_status_counts <- list(analysed = 2)
fails(contrast_rows(lead, comparisons, frequencies, bad_summary), "Mismatched status totals accepted")

folder <- tempfile("R_contrast_test_"); dir.create(folder)
draw_haplotype_contrasts(rows, summary, file.path(folder, "synthetic_tested.pdf"), "pdf")
draw_regional_haplotype_contrasts(rows, summary, file.path(folder, "synthetic_regional.pdf"), "pdf")
draw_haplotype_contrasts(contrast_rows(leads, comparisons[0, ], frequencies[0, ], summary_for(leads)), summary_for(leads), file.path(folder, "synthetic_unavailable.pdf"), "pdf")
assert(all(file.info(list.files(folder, full.names = TRUE))$size > 1000), "Contrast graphic did not render")
assert(any(grepl("descriptive-only differences", contrast_interpretation(rows, summary))), "Interpretation lost descriptive differences")

# Export requires a completed, stable source and fingerprints every input, including completion.
run <- file.path(folder, "source"); dir.create(run)
write.table(lead, file.path(run, "lead_loci.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(comparisons, file.path(run, "haplotype_comparisons.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(frequencies, file.path(run, "haplotype_group_frequencies.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
jsonlite::write_json(summary, file.path(run, "run_summary.json"), auto_unbox = TRUE)
fails(export_haplotype_contrasts(run, file.path(folder, "incomplete")), "Incomplete analysis exported")
writeLines("synthetic", file.path(run, "COMPLETED.txt"))
unstable <- summary; unstable$inputs_stable_during_read <- FALSE
jsonlite::write_json(unstable, file.path(run, "run_summary.json"), auto_unbox = TRUE)
fails(export_haplotype_contrasts(run, file.path(folder, "unstable")), "Unstable source exported")
jsonlite::write_json(summary, file.path(run, "run_summary.json"), auto_unbox = TRUE)
out <- file.path(folder, "output"); export_haplotype_contrasts(run, out)
manifest <- jsonlite::fromJSON(file.path(out, "manifest.json"))
assert(identical(manifest$source_run_path, normalizePath(run)) && any(grepl(normalizePath(run), readLines(file.path(out, "index.md")), fixed = TRUE)), "Full source run path missing")
assert(length(manifest$source_sha256) == 5L && "COMPLETED.txt" %in% names(manifest$source_sha256), "Completion fingerprint missing")
assert(manifest$comparison_status_counts$tested == 1 && manifest$comparison_status_counts$descriptive_only == 1, "Manifest statuses disagree")
assert(file.exists(file.path(out, "regional", "L0001", "haplotype_contrasts.png")), "Regional picture missing")
fails(export_haplotype_contrasts(run, out), "Existing output was overwritten")
# Simulate a source changing during rendering: no success manifest may be written.
original_draw <- draw_haplotype_contrasts
changed <- FALSE
draw_haplotype_contrasts <- function(...) {
  original_draw(...)
  if (!changed) {
    cat("changed\n", file = file.path(run, "COMPLETED.txt"), append = TRUE)
    changed <<- TRUE
  }
}
changed_out <- file.path(folder, "changed_source")
fails(export_haplotype_contrasts(run, changed_out), "Source changed during rendering was accepted")
assert(!file.exists(file.path(changed_out, "manifest.json")), "Changed source received a success manifest")
draw_haplotype_contrasts <- original_draw
unlink(folder, recursive = TRUE)
cat("PASS: native R observed carrier contrasts, inference eligibility, orientation, zero/missing distinction, schema validation, source fingerprints and overview/regional rendering.\n")
