#!/usr/bin/env Rscript
entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(entry)), "finemap_external.R"))
assert <- function(value, message) if (!isTRUE(value)) stop(message, call. = FALSE)
fails <- function(expression, message, pattern = NULL) {
  result <- tryCatch(force(expression), error = identity)
  assert(inherits(result, "error"), message)
  if (!is.null(pattern)) assert(grepl(pattern, conditionMessage(result), fixed = TRUE), paste(message, conditionMessage(result)))
}
folder <- tempfile("external_synthetic_"); dir.create(folder)
on.exit_cleanup <- function() unlink(folder, recursive = TRUE)
input <- file.path(folder, "input"); dir.create(input); dir.create(file.path(input, "SYNTHETIC"))
R <- toeplitz(.6^(0:7)); z <- R[,2]*8+R[,5]*.3
variants <- data.frame(variant_id = paste0("synthetic_v", seq_len(8)), chrom = "chr1", pos = seq(1100, 9500, length.out = 8),
  ref = c("A", "A", "AT", "G", "C", "G", "A", "C"), alt = c("G", "AT", "A", "A", "T", "T", "C", "G"),
  variant_type = c("SNP", "INDEL", "INDEL", rep("SNP", 5)), beta = z*.01, se = .01,
  pval = 2*pnorm(-abs(z)), af = .3, gwas_n = 30000, reference_af = .32)
source_entry <- list(locus_id = "SYNTHETIC", variants_file = "SYNTHETIC/matched_variants.tsv", ld_file = "SYNTHETIC/ld.tsv.gz",
  window_start = 1000, window_end = 11000, ld_variant_ids = as.list(variants$variant_id))
manifest <- list(synthetic = TRUE, n_cases = 10000, n_controls = 20000,
  provenance = list(gwas_source = "SYNTHETIC regression summaries", ld_source = "SYNTHETIC Toeplitz correlation", genome_build = "GRCh38",
    effect_allele = "ALT", ld_dosage_allele = "ALT", allele_alignment = "exact_normalized_REF_ALT", reference_missingness = "complete_calls_only"),
  loci = list(source_entry))
write_fixture <- function(v = variants, r = R, entry = source_entry, m = manifest) {
  external_write_table(v, file.path(input, entry$variants_file))
  con <- gzfile(file.path(input, entry$ld_file), "wt"); write.table(r, con, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, na = "NA"); close(con)
  entry$variants_sha256 <- external_sha256(file.path(input, entry$variants_file)); entry$ld_sha256 <- external_sha256(file.path(input, entry$ld_file))
  m$loci <- list(entry); external_json(m, file.path(input, "input_manifest.json"))
  external_manifest(input)
}
load_fixture <- function(m) external_load_locus(input, m$data$loci[[1]], m$total_n)
m <- write_fixture(); good <- load_fixture(m)
assert(identical(good$variants$variant_type, c("SNP", "insertion", "deletion", rep("SNP", 5))), "SNP and indel types were lost")
assert(max(abs(good$R-R)) < 1e-12, "LD was changed or shrunk")
assert(abs(m$effective_n-4/(1/10000+1/20000)) < 1e-9, "Effective N formula is wrong")
complex <- variants; complex$ref[2] <- "AT"; complex$alt[2] <- "ACG"
assert(load_fixture(write_fixture(complex))$variants$variant_type[2] == "complex_indel", "Complex length-changing allele mislabeled pure insertion")
complex$ref[2] <- "A"; complex$alt[2] <- "TG"
assert(load_fixture(write_fixture(complex))$variants$variant_type[2] == "complex_indel", "Unanchored replacement mislabeled pure insertion")
rich <- manifest; rich$provenance$gwas_source <- list(url = "https://example.invalid/SYNTHETIC", synthetic = TRUE)
rich$provenance$ld_source <- list(reference_samples = 503, synthetic = TRUE)
assert(is.list(write_fixture(m = rich)$data$provenance$gwas_source), "Structured source provenance rejected")
m <- write_fixture()
# Aligned known indel signal is fitted jointly with SNPs and yields a credible set.
out <- file.path(folder, "successful")
run <- run_external_finemapping(input, out, workers = 1, locus = "SYNTHETIC")
pip <- read.delim(file.path(out, "SYNTHETIC", "variant_pip.tsv"), na.strings = "NA")
summary <- jsonlite::fromJSON(file.path(out, "SYNTHETIC", "summary.json"))
cs <- read.delim(file.path(out, "SYNTHETIC", "credible_sets.tsv"))
assert(summary$status == "exploratory_window_limited" && summary$converged && summary$n_insertions == 1 && summary$n_deletions == 1, "Synthetic fit or type counts failed")
assert(pip$pip[2] > .95 && which.max(pip$pip) == 2 && pip$variant_type[2] == "insertion", "Known synthetic indel signal was not recovered")
assert(summary$sensitivity_status == "converged" && is.finite(summary$sensitivity_max_abs_pip_difference), "Effective-N sensitivity missing")
assert(abs(summary$primary_n-m$effective_n) < 1e-8 && summary$sensitivity_n == m$total_n && run$primary_n == m$effective_n, "Primary effective-N and sensitivity total-N roles swapped")
assert(nrow(cs) >= 1 && any(cs$boundary_sensitive) && any(grepl("synthetic_v2", cs$variant_ids, fixed = TRUE)), "Boundary credible set or member mapping wrong")
assert(!summary$validated_causal_candidates && run$analysis == "exploratory_window_limited_external_SuSiE_RSS" && !run$cohort_haplotype_inference_modified, "Exploratory external scope lost")
assert(run$script_sha256 == external_sha256(file.path(out, "finemap_external.R")), "Exact fitting script provenance missing")
assert(all(file.exists(file.path(out, "SYNTHETIC", c("model_total_n.rds", "model_effective_n.rds", "external_finemapping.png", "diagnostics.tsv")))), "Required model artifacts missing")
assert(file.info(file.path(out, "SYNTHETIC", "external_finemapping.png"))$size > 1000, "Plot is empty")
fails(run_external_finemapping(input, out, 1), "Output preservation failed", "new or empty")
fails(run_external_finemapping(input, file.path(folder, "unknown"), 1, "ABSENT"), "Unknown pilot locus accepted")
# PIPs are invariant to a consistent joint permutation, and alignment errors fail.
order <- c(8, 2, 6, 3, 7, 1, 4, 5)
permuted <- source_entry; permuted$ld_variant_ids <- as.list(variants$variant_id[order])
pm <- write_fixture(variants[order, ], R[order, order], permuted)
pg <- load_fixture(pm)
pfit <- external_fit(pg$variants$beta/pg$variants$se, pg$R, pm$effective_n)
ppip <- external_usable_fit(pfit, 8)
assert(max(abs(ppip[match(variants$variant_id, pg$variants$variant_id)]-pip$pip)) < 1e-7, "Joint SNP/indel ordering changed model probabilities")
wrong_order <- source_entry; wrong_order$ld_variant_ids <- rev(wrong_order$ld_variant_ids)
fails(load_fixture(write_fixture(entry = wrong_order)), "Mismatched unlabeled LD order accepted", "LD variant order")
# Input filtering and typing are verified again rather than trusted implicitly.
for (field in c("beta", "se", "af", "reference_af", "pos", "gwas_n")) {
  bad <- variants; bad[[field]][1] <- NA
  fails(load_fixture(write_fixture(bad)), paste("Missing", field, "accepted"))
}
bad <- variants; bad$variant_id[2] <- bad$variant_id[1]
fails(load_fixture(write_fixture(bad)), "Duplicate IDs accepted")
bad <- variants; bad$variant_type[2] <- "SNP"
fails(load_fixture(write_fixture(bad)), "Indel mislabeled SNP accepted", "Variant type disagrees")
bad <- variants; bad$se[1] <- 0
fails(load_fixture(write_fixture(bad)), "Zero SE accepted")
bad <- variants; bad$af[1] <- .001
fails(load_fixture(write_fixture(bad)), "Rare variant bypassed common-variant filter")
bad <- variants; bad$reference_af[1] <- .45
fails(load_fixture(write_fixture(bad)), "Allele frequency mismatch accepted")
bad <- variants; bad$gwas_n[1] <- 29999
fails(load_fixture(write_fixture(bad)), "Inconsistent nominal study N accepted")
bad <- manifest; bad$provenance$ld_dosage_allele <- "REF"
fails(write_fixture(m = bad), "REF-coded LD provenance accepted")
m <- write_fixture(); cat("tamper\n", file = file.path(input, source_entry$variants_file), append = TRUE)
fails(load_fixture(m), "Changed source hash accepted", "SHA256 mismatch")
for (kind in c("missing", "asymmetric", "diagonal", "nonpsd", "dimension")) {
  bad <- R
  if (kind == "missing") bad[1, 2] <- NA
  if (kind == "asymmetric") bad[1, 2] <- .1
  if (kind == "diagonal") bad[1, 1] <- .9
  if (kind == "nonpsd") { bad <- matrix(-.8, 8, 8); diag(bad) <- 1 }
  if (kind == "dimension") bad <- bad[-1, -1]
  fails(load_fixture(write_fixture(r = bad)), paste("Invalid LD accepted:", kind))
}
# Valid rank deficiency from shared LD is accepted, not repaired or shrunk.
singular <- R; singular[8, ] <- singular[7, ]; singular[, 8] <- singular[, 7]
sg <- load_fixture(write_fixture(r = singular))
assert(abs(sg$minimum_eigenvalue) < 1e-7 && sg$R[7, 8] == 1, "Valid semidefinite LD rejected or changed")
# Severe incompatibility produces diagnostics and external estimates, never PIPs.
bad <- variants; bad$beta <- rep(c(.08, -.08), 4)
bad$pval <- 2*pnorm(-abs(bad$beta/bad$se)); m <- write_fixture(bad)
blocked_out <- file.path(folder, "blocked")
run_external_finemapping(input, blocked_out, 1)
blocked <- jsonlite::fromJSON(file.path(blocked_out, "SYNTHETIC", "summary.json"))
blocked_pip <- read.delim(file.path(blocked_out, "SYNTHETIC", "variant_pip.tsv"), na.strings = "NA")
assert(blocked$status == "unavailable_severe_LD_incompatibility" && blocked$s >= .5, "Predeclared severe mismatch block missing")
assert(!any(grepl("jointly fitted", readLines(file.path(blocked_out, "SYNTHETIC", "INTERPRETATION.txt")), fixed = TRUE)), "Blocked interpretation incorrectly claims a fit")
assert(all(is.na(blocked_pip$pip)) && all(is.na(blocked_pip$pip_total_n)) && !file.exists(file.path(blocked_out, "SYNTHETIC", "model_total_n.rds")), "Blocked fit leaked PIPs")
assert(nrow(read.delim(file.path(blocked_out, "SYNTHETIC", "credible_sets.tsv"))) == 0, "Blocked fit leaked credible sets")
# An unconverged fit must not leak its available internal numeric probabilities.
real_fit <- external_fit
external_fit <- function(...) { result <- real_fit(...); result$converged <- FALSE; result }
m <- write_fixture(); nonconverged_out <- file.path(folder, "nonconverged")
run_external_finemapping(input, nonconverged_out, 1)
failed <- jsonlite::fromJSON(file.path(nonconverged_out, "SYNTHETIC", "summary.json"))
assert(failed$status == "unavailable_fit" && !failed$converged, "Nonconvergence status lost")
assert(all(is.na(read.delim(file.path(nonconverged_out, "SYNTHETIC", "variant_pip.tsv"))$pip)), "Unconverged PIPs were retained")
external_fit <- real_fit
# Retain finite converged diagnostic-warning fits only under an explicit label.
real_diagnostics <- external_diagnostics
external_diagnostics <- function(...) { result <- real_diagnostics(...); result$s <- .2; result$warning <- TRUE; result }
warning_out <- file.path(folder, "warning"); run_external_finemapping(input, warning_out, 1)
warning_summary <- jsonlite::fromJSON(file.path(warning_out, "SYNTHETIC", "summary.json"))
assert(warning_summary$status == "exploratory_with_LD_warning" && !warning_summary$validated_causal_candidates, "Warning fit mislabeled as validated")
external_diagnostics <- real_diagnostics
# File mutation after read discards both primary and sensitivity results.
external_fit <- function(...) {
 result <- real_fit(...)
 cat("changed\n", file = file.path(input, source_entry$variants_file), append = TRUE)
 result
}
changed_out <- file.path(folder, "changed"); run_external_finemapping(input, changed_out, 1)
changed <- jsonlite::fromJSON(file.path(changed_out, "SYNTHETIC", "summary.json"))
assert(changed$status == "unavailable_input_changed" && !file.exists(file.path(changed_out, "SYNTHETIC", "model_total_n.rds")), "Changed input retained model")
assert(all(is.na(read.delim(file.path(changed_out, "SYNTHETIC", "variant_pip.tsv"))$pip)), "Changed input retained probabilities")
external_fit <- real_fit
# Invalid inputs still get explicit locus status and an empty PIP table.
m <- write_fixture(); cat("changed\n", file = file.path(input, source_entry$variants_file), append = TRUE)
invalid_out <- file.path(folder, "invalid"); run_external_finemapping(input, invalid_out, 1)
invalid <- jsonlite::fromJSON(file.path(invalid_out, "SYNTHETIC", "summary.json"))
assert(invalid$status == "unavailable_input_validation" && nrow(read.delim(file.path(invalid_out, "SYNTHETIC", "variant_pip.tsv"))) == 0, "Invalid input manufactured estimates")
# Independent batch workers preserve row order and reproduce the pilot model.
m <- write_fixture(); twin <- m$data$loci[[1]]; twin$locus_id <- "SYNTHETIC_B"
m$data$loci[[2]] <- twin; external_json(m$data, file.path(input, "input_manifest.json"))
parallel_out <- file.path(folder, "parallel"); batch <- run_external_finemapping(input, parallel_out, 2)
parallel_pip <- read.delim(file.path(parallel_out, "SYNTHETIC_B", "variant_pip.tsv"))
assert(batch$loci_attempted == 2 && max(abs(parallel_pip$pip-pip$pip)) < 1e-12, "Parallel batch disagrees with pilot")
assert(identical(read.delim(file.path(parallel_out, "locus_summary.tsv"))$locus_id, c("SYNTHETIC", "SYNTHETIC_B")), "Batch lost manifest locus order")
on.exit_cleanup()
cat("PASS: synthetic joint SNP/indel SuSiE signal, LD and hash/order validation, effective-N sensitivity, boundary CS, mismatch/nonconvergence/input-change gates, private plots and preserved runs.\n")
