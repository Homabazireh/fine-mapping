#!/usr/bin/env Rscript
entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(entry)), "compare_haplotypes.R"))
source(file.path(dirname(normalizePath(entry)), "plot_haplotype_contrasts.R"))
assert <- function(x, message) if (!isTRUE(x)) stop(message, call. = FALSE)
near <- function(a, b) assert(isTRUE(all.equal(a, b, tolerance = 1e-12)), "Numeric contrast mismatch")
leads <- data.frame(locus_id = c("L0001", "L0002"), chrom = c("chr1", "chr2"), lead_pos = c("100", "200"),
                    status = c("lead_below_min_call_rate", "palindromic_lead_orientation_unconfirmed"))
empty <- function(columns) as.data.frame(setNames(rep(list(character()), length(columns)), columns), stringsAsFactors = FALSE)
comparisons <- empty(COMP_FIELDS); frequencies <- empty(FREQ_FIELDS)
summary <- list(group_counts = list(T2D = 20, Prediabetes = 24, Control = 30), selected_loci = 2,
                synthetic = TRUE, parameters = list(min_call_rate = .9))
rows <- contrast_rows(leads, comparisons, frequencies, summary)
assert(length(rows) == 4L, "Skipped loci must appear in both contrasts")
assert(all(vapply(rows, function(r) r$comparison_status == "not_estimable" && is.na(r$carrier_frequency_difference) && is.na(r$CI95_low) && is.na(r$CI95_high), logical(1))), "Unavailable contrasts became zero effects")
assert(identical(vapply(rows, `[[`, character(1), "reason"), rep(leads$status, each = 2)), "Exclusion reasons changed")

lead <- leads[1, , drop = FALSE]; lead$status <- "analysed"
source_rows <- list(list(locus_id = "L0001", haplotype_id = "L0001_H001", group_a = "Control", group_b = "T2D",
  n_resolved_a = "30", n_resolved_b = "20", carriers_a = "6", carriers_b = "8", carrier_frequency_a = "0.2", carrier_frequency_b = "0.4",
  test_status = "tested_unadjusted", carrier_frequency_difference_a_minus_b = "-0.2", carrier_difference_CI95_low = "-0.38",
  carrier_difference_CI95_high = "-0.03", p_fisher_two_sided = "0.02", q_BH_all_pairwise_tests = "0.08"),
  list(locus_id = "L0001", haplotype_id = "L0001_H001", group_a = "Prediabetes", group_b = "Control",
  n_resolved_a = "24", n_resolved_b = "30", carriers_a = "6", carriers_b = "6", carrier_frequency_a = "0.25", carrier_frequency_b = "0.2",
  test_status = "not_tested_reference_provenance", carrier_frequency_difference_a_minus_b = "NA", carrier_difference_CI95_low = "NA",
  carrier_difference_CI95_high = "NA", p_fisher_two_sided = "NA", q_BH_all_pairwise_tests = "NA"))
comparisons <- do.call(rbind, lapply(source_rows, as.data.frame, stringsAsFactors = FALSE))
frequencies <- data.frame(locus_id = "L0001", haplotype_id = "L0001_H001", group = c("T2D", "Prediabetes", "Control"),
                          haplotype_frequency = c("0.3", "0.2", "0.1"))
rows <- contrast_rows(lead, comparisons, frequencies, summary)
disease <- rows[[1]]; prediabetes <- rows[[2]]
near(disease$carrier_frequency_difference, .2)
near(disease$CI95_low, .03); near(disease$CI95_high, .38)
near(disease$q_BH_all_pairwise_tests, .08)
near(disease$n_resolved_case, 20); near(disease$carriers_case, 8)
near(disease$haplotype_frequency_case_descriptive, .3)
assert(prediabetes$comparison_status == "not_estimable" && is.na(prediabetes$carrier_frequency_difference), "Provenance gate bypassed")
assert(prediabetes$reason == "not_tested_reference_provenance", "Provenance reason missing")
forward <- comparisons
forward[1, c("group_a", "group_b")] <- c("T2D", "Control")
forward[1, c("carrier_frequency_difference_a_minus_b", "carrier_difference_CI95_low", "carrier_difference_CI95_high")] <- c("0.2", "0.03", "0.38")
near(contrast_rows(lead, forward, frequencies, summary)[[1]]$carrier_frequency_difference, disease$carrier_frequency_difference)
bad <- comparisons; bad$carrier_difference_CI95_high[1] <- "NA"
assert(inherits(try(contrast_rows(lead, bad, frequencies, summary), silent = TRUE), "try-error"), "Tested contrast with missing uncertainty accepted")
missing <- contrast_rows(lead, comparisons[0, ], frequencies, summary)
assert(all(vapply(missing, function(r) r$reason == "pairwise_result_unavailable" && is.na(r$carrier_frequency_difference), logical(1))), "Missing pair became an effect")
no_haps <- contrast_rows(lead, comparisons[0, ], frequencies[0, ], summary)
assert(all(vapply(no_haps, function(r) r$reason == "no_resolved_haplotypes", logical(1))), "Empty resolved panel reason lost")
folder <- tempfile("R_contrast_test_"); dir.create(folder)
draw_haplotype_contrasts(rows, summary, file.path(folder, "synthetic_tested.pdf"), "pdf")
draw_haplotype_contrasts(contrast_rows(leads, comparisons[0, ], frequencies[0, ], summary), summary, file.path(folder, "synthetic_unavailable.pdf"), "pdf")
assert(all(file.info(list.files(folder, full.names = TRUE))$size > 1000), "Contrast graphic did not render")
unlink(folder, recursive = TRUE)
cat("PASS: native R contrast orientation, reversed CIs, preserved q, missingness, provenance and both graphic branches.\n")
