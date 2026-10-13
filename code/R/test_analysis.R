#!/usr/bin/env Rscript
# Tests generate software fixtures in R; no Python execution or private input needed.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script)), "compare_haplotypes.R"))

assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
near <- function(a, b, tolerance = 1e-10) assert(isTRUE(all.equal(unname(a), unname(b), tolerance = tolerance)), paste("Numeric mismatch:", paste(a, collapse = ","), "vs", paste(b, collapse = ",")))
fails <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  assert(inherits(error, "error") && grepl(pattern, conditionMessage(error)), paste("Expected error:", pattern))
}

fmt <- c("GT", "DP", "GQ", "PS")
call <- function(s) parse_call(s, fmt, 2, 10, 20)
assert(is.null(call("./.:20:40:1")$gt), "Missing GT was filled")
assert(is.null(call("0/.:20:40:1")$gt), "Partial GT was filled")
assert(identical(call("0/0:20:40:.")$gt, c(0L, 0L)), "Reference GT lost")
assert(call("0|1:9:40:1")$reason == "low_or_invalid_DP", "DP threshold")
assert(call("0|1:20:19:1")$reason == "low_or_invalid_GQ", "GQ threshold")
assert(call("0|2:20:40:1")$reason == "invalid_allele_index", "Invalid allele index")
assert(call("0|1|:20:40:1")$reason == "not_diploid", "Malformed trailing ploidy accepted")
assert(parse_call("0/1:FAIL", c("GT", "FT"), 2, 0, 0)$reason == "genotype_filter_failed", "FT filtering")
assert(identical(parse_call("0/1", "GT", 2, 0, 0)$gt, c(0L, 1L)), "Disabled DP/GQ require missing fields")
assert(is.null(classify_variant("A", "<DEL>", 49)), "Symbolic variant included")
assert(classify_variant("A", c("C", "AT"), 49) == "SNP/indel", "Multiallelic site broken apart")
assert(resolve_diplotype(list(call("0/1:20:40:."), call("1/1:20:40:.")))$status == "single_het_resolved", "One heterozygote not resolved")
assert(resolve_diplotype(list(call("0/1:20:40:."), call("0/1:20:40:.")))$status == "unphased_heterozygote", "Unphased sites stitched")
assert(resolve_diplotype(list(call("0|1:20:40:1"), call("0|1:20:40:2")))$status == "different_phase_sets", "Phase boundary stitched")
assert(resolve_diplotype(list(call("0|1:20:40:."), call("0|1:20:40:.")))$status == "missing_PS", "Missing phase set inferred")
forward <- resolve_diplotype(list(call("0|1:20:40:1"), call("1|0:20:40:1")))$pair
reverse <- resolve_diplotype(list(call("1|0:20:40:1"), call("0|1:20:40:1")))$pair
assert(identical(forward, reverse), "Global phase flip changes diplotype")
near(decimal_log10(c("0", "1e-500", "2.5E-400", "0.001", ".2", "1")), c(-Inf, -500, log10(2.5)-400, -3, log10(.2), 0))
fails(decimal_log10("1.1"), "outside")
fails(decimal_log10("1.0000000000000001"), "outside")
assert(identical(decimal_at_most(c("1e-8", "1.0000000000000001e-8", "9.9999999999999999e-9"), "1e-8"), c(TRUE, FALSE, TRUE)), "Exact p-value threshold precision")
test <- carrier_test(1, 10, 11, 14)
near(test$p_fisher_two_sided, 0.0027594561852200836)
near(test$conditional_odds_ratio_a_vs_b, 0.037209084832381056)
near(carrier_test(1, 3, 1, 4)$CI95_high, 156.4984076271612, 1e-10)
assert(test$carrier_difference_CI95_low < test$carrier_frequency_difference_a_minus_b &&
       test$carrier_difference_CI95_high > test$carrier_frequency_difference_a_minus_b, "Difference CI excludes estimate")
overall <- overall_test(c(3, 8, 6), c(10, 10, 10))
near(overall$p, 0.10299286447002098)
near(overall$p, fisher.test(cbind(c(3, 8, 6), c(7, 2, 4)))$p.value)
near(overall_test(c(1, 11), c(10, 14))$p, test$p_fisher_two_sided)
near(overall_test(c(6, 3, 8), c(10, 10, 10))$p, overall$p)

folder <- tempfile("native_R_test_"); dir.create(folder)
inputs <- file.path(folder, "inputs"); dir.create(inputs)
samples <- paste0(rep(c("D", "P", "C"), each = 4), rep(1:4, 3))
groups <- rep(c("T2D", "Prediabetes", "Control"), each = 4)
write.table(data.frame(sample = samples, group = groups), file.path(inputs, "samples.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
writeLines(c("Chromsome\tPosition\tEffectAllele\tNonEffectAllele\tBeta\tSE\tPval\tliftover_strand",
             "1\t10000\tG\tA\t0.2\t0.03\t1e-500\t+",
             "1\t200010\tT\tC\t0.2\t0.04\t1e-10\t+"), file.path(inputs, "gwas.tsv"))
header <- c("##fileformat=VCFv4.2", "##reference=GRCh38", "##contig=<ID=chr1,length=1000000>",
            paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", samples), collapse = "\t"))
haps <- rbind(c(0,0,0,0), c(1,1,1,1), c(1,0,1,0))
pairs <- rbind(c(2,2), c(2,1), c(2,3), c(2,1), c(2,1), c(2,3), c(1,1), c(1,3), c(1,1), c(1,1), c(1,3), c(1,1))
records <- character()
for (j in 1:4) {
  pos <- c(9990, 9995, 10000, 10005)[j]; ref <- c("C", "AT", "A", "G")[j]; alt <- c("T", "A", "G", "C")[j]
  cells <- vapply(seq_along(samples), function(i) {
    ps <- 1000+i + if (i == 4 && j == 4) 100 else 0
    if (i == 8 && j == 2) "./.:.:.:." else paste0(haps[pairs[i,1],j], "|", haps[pairs[i,2],j], ":25:40:", ps)
  }, character(1))
  records <- c(records, paste(c("chr1", pos, ".", ref, alt, "60", "PASS", ".", "GT:DP:GQ:PS", cells), collapse = "\t"))
}
for (j in 1:3) {
  cells <- vapply(0:11, function(i) paste0(i %% 2, "|", (i %/% 2) %% 2, ":25:40:", 2000+i), character(1))
  records <- c(records, paste(c("chr1", c(200005,200010,200018)[j], ".", c("T","C","G")[j], c("TA","T","A")[j], "60", "PASS", ".", "GT:DP:GQ:PS", cells), collapse = "\t"))
}
writeLines(c(header, records), file.path(inputs, "cohort.vcf"))
args <- defaults(); args$vcf <- file.path(inputs, "cohort.vcf"); args$gwas <- file.path(inputs, "gwas.tsv")
args$samples <- file.path(inputs, "samples.tsv"); args$out <- file.path(folder, "requested")
args$genome_build <- "GRCh38"; args$synthetic <- TRUE; args$figure <- FALSE; args$test_associations <- TRUE
args$min_test_samples <- 2; args$min_resolved_fraction <- .7; args$min_carriers <- 1
summary <- run_analysis(args)
assert(summary$haplotype_carrier_tests_performed == 0, "Unconfirmed provenance allowed tests")
pairwise <- read_tsv(file.path(args$out, "haplotype_comparisons.tsv"))
assert(all(pairwise$test_status == "not_tested_reference_provenance"), "Provenance reason lost")
assert(all(pairwise$p_fisher_two_sided == "NA"), "Unconfirmed provenance has p-values")
assignment <- read_tsv(file.path(args$out, "sample_haplotypes.tsv"))
assert(assignment$status[assignment$locus_id == "L0001" & assignment$sample == "D4"] == "different_phase_sets", "Synthetic PS boundary")
assert(assignment$status[assignment$locus_id == "L0001" & assignment$sample == "P4"] == "missing_or_low_quality", "Synthetic missing call")
freq <- read_tsv(file.path(args$out, "haplotype_group_frequencies.tsv"))
for (key in unique(paste(freq$locus_id, freq$group))) {
  rows <- freq[paste(freq$locus_id, freq$group) == key, ]
  near(sum(as.numeric(rows$haplotype_frequency)), 1)
  near(sum(as.numeric(rows$haplotype_copies)), unique(as.numeric(rows$resolved_chromosomes)))
  assert(all(as.numeric(rows$carriers) <= as.numeric(rows$n_resolved)), "Carrier counts exceed people")
}
fails(run_analysis(args), "not empty")
args$out <- file.path(folder, "confirmed"); args$reference_genotypes_confirmed <- TRUE
summary <- run_analysis(args)
assert(summary$haplotype_carrier_tests_performed > 0 && summary$haplotype_overall_tests_performed > 0, "Confirmed synthetic association tests missing")
pw <- read_tsv(file.path(args$out, "haplotype_comparisons.tsv")); ov <- read_tsv(file.path(args$out, "haplotype_overall_comparisons.tsv"))
for (spec in list(list(x = pw, p = "p_fisher_two_sided", q = "q_BH_all_pairwise_tests"), list(x = ov, p = "p_fisher_overall", q = "q_BH_all_overall_tests"))) {
  valid <- spec$x$test_status == "tested_unadjusted"
  near(as.numeric(spec$x[[spec$q]][valid]), p.adjust(as.numeric(spec$x[[spec$p]][valid]), "BH"))
}
args$out <- file.path(folder, "no_panel"); args$min_call_rate <- 1; args$min_dp <- 100
summary <- run_analysis(args)
assert(nrow(read_tsv(file.path(args$out, "haplotype_group_frequencies.tsv"))) == 0, "Invalid lead yielded frequencies")
assert(nrow(read_tsv(file.path(args$out, "haplotype_overall_comparisons.tsv"))) == 0, "No panel yielded tests")
unlink(folder, recursive = TRUE)
cat("PASS: native R genotype, phase, exact-test, interval, correction, provenance and end-to-end tests.\n")
