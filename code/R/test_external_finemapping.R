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
  invisible(external_manifest(input))
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
# mixsqp/irlba uses randomized truncated SVD inside kriging: seed locally and
# preserve the caller's RNG kind/state on success, failure and uninitialized use.
test_diagnostic_rng <- function() {
  old_kind <- RNGkind(); had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({ do.call(RNGkind, as.list(old_kind)); if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv) })
  dR <- toeplitz(.7^(0:63)); dz <- dR[, 8]*8+dR[, 25]*3
  dv <- data.frame(variant_id = paste0("diagnostic_synthetic_", 1:64), variant_type = "SNP", beta = dz*.01, se = .01)
  set.seed(123L, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
  caller_kind <- RNGkind(); caller_seed <- .Random.seed
  first <- external_diagnostics(dv, dR, 30000)
  assert(identical(RNGkind(), caller_kind) && identical(.Random.seed, caller_seed), "Successful diagnostic changed caller RNG kind/state")
  set.seed(987L, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  other_seed <- .Random.seed; second <- external_diagnostics(dv, dR, 30000)
  assert(identical(first, second) && identical(.Random.seed, other_seed), "Diagnostics vary with caller seed or RNG algorithm")
  expected_rng <- external_with_diagnostic_rng(c(runif(3), rnorm(3)))
  set.seed(456L); assert(identical(expected_rng, external_with_diagnostic_rng(c(runif(3), rnorm(3)))), "Local RNG helper failed to reset its declared seed")
  caller_kind <- RNGkind(); caller_seed <- .Random.seed
  fails(external_with_diagnostic_rng({ runif(2); stop("synthetic diagnostic error") }), "Diagnostic error was swallowed", "synthetic diagnostic error")
  assert(identical(RNGkind(), caller_kind) && identical(.Random.seed, caller_seed), "Failed diagnostic changed caller RNG kind/state")
  rm(".Random.seed", envir = .GlobalEnv)
  absent <- external_diagnostics(dv, dR, 30000)
  assert(identical(first, absent) && !exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE) && identical(RNGkind(), caller_kind), "Diagnostic created a persistent RNG seed when none existed")
  fails(external_with_diagnostic_rng({ rnorm(2); stop("synthetic absent-seed error") }), "Absent-seed error was swallowed")
  assert(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE) && identical(RNGkind(), caller_kind), "Diagnostic error created a persistent RNG seed")
  parallel_diagnostics <- parallel::mclapply(c(11L, 912L), function(seed) { set.seed(seed); external_diagnostics(dv, dR, 30000) }, mc.cores = 2L)
  assert(all(vapply(parallel_diagnostics, identical, logical(1), first)), "Diagnostic output depends on worker RNG state")
}
test_diagnostic_rng()
for (name in c("diagnostics.tsv", "diagnostics_total_n.tsv")) {
  assert(identical(read.delim(file.path(out, "SYNTHETIC", name)), read.delim(file.path(parallel_out, "SYNTHETIC_B", name))), "Actual single-worker and parallel-worker diagnostic tables disagree")
}
assert(identical(run$diagnostic_rng, EXTERNAL_DIAGNOSTIC_RNG) && all(c("mixsqp", "irlba") %in% names(run$packages)), "Diagnostic RNG policy or numerical dependency provenance missing")

# Explicit JSON defaults reproduce the same model, with exact configuration provenance.
m <- write_fixture()
config_dir <- file.path(folder, "configuration files"); dir.create(config_dir)
config_path <- file.path(config_dir, "model options.json")
config <- list(paths = list(input_dir = "../input", out = "../configured"),
  analysis = list(workers = 1, locus = "SYNTHETIC", max_effects = 10, coverage = .95, max_iter = 200, variant_types = "both"),
  plot = list(width = 1800, height = 1150, dpi = 150, formats = list("png")))
external_json(config, config_path)
configured <- run_external_finemapping_config(config_path)
configured_out <- file.path(folder, "configured")
configured_pip <- read.delim(file.path(configured_out, "SYNTHETIC", "variant_pip.tsv"))
assert(max(abs(configured_pip$pip-pip$pip)) < 1e-12 && max(abs(configured_pip$pip_total_n-pip$pip_total_n)) < 1e-12, "Explicit default configuration changed PIPs")
assert(identical(read.delim(file.path(configured_out, "SYNTHETIC", "credible_sets.tsv")), cs), "Explicit default configuration changed credible sets")
assert(configured$configuration$sha256 == external_sha256(config_path) && external_sha256(file.path(configured_out, configured$configuration$snapshot)) == external_sha256(config_path), "Exact configuration snapshot/hash lost")
assert(configured$effective_settings$paths$input_dir == normalizePath(input) && configured$effective_settings$paths$out == normalizePath(configured_out), "JSON paths did not resolve from its directory")
assert(identical(configured$thresholds, run$thresholds), "Default configuration altered fixed QC thresholds")
assert(identical(external_resolve_options(config = file.path(configured_out, "effective_settings.json"))$options$plot_formats, "png"), "Effective settings are not reusable configuration")

# Actual CLI overrides reach both models, credible sets and raster/vector output.
custom_out <- file.path(folder, "custom")
custom <- external_cli(c("--config", config_path, "--out", custom_out, "--workers", "1", "--max-effects", "2", "--coverage", "0.8",
  "--max-iter", "25", "--plot-width", "1200", "--plot-height", "800", "--plot-dpi", "100", "--plot-formats", "png,pdf,svg"))
custom_summary <- jsonlite::fromJSON(file.path(custom_out, "SYNTHETIC", "summary.json"))
assert(custom$max_effects == 2 && custom$coverage == .8 && custom$max_iter == 25 && custom_summary$L == 2 && custom_summary$max_iter == 25, "CLI overrides did not reach output settings")
for (name in c("model_effective_n.rds", "model_total_n.rds")) {
  fit <- readRDS(file.path(custom_out, "SYNTHETIC", name))
  assert(nrow(fit$alpha) == 2 && fit$sets$requested_coverage == .8 && fit$niter <= 25, "Model did not use requested effect count, coverage or iterations")
}
png_dimensions <- function(path) {
  con <- file(path, "rb"); on.exit(close(con)); header <- readBin(con, "raw", 24L)
  c(sum(as.integer(header[17:20])*256^(3:0)), sum(as.integer(header[21:24])*256^(3:0)))
}
assert(identical(png_dimensions(file.path(custom_out, "SYNTHETIC", "external_finemapping.png")), c(1200, 800)), "PNG pixel dimensions ignored")
assert(identical(png_dimensions(file.path(configured_out, "SYNTHETIC", "external_finemapping.png")), c(1800, 1150)), "Default PNG dimensions changed")
vector_pdf <- suppressWarnings(readLines(file.path(custom_out, "SYNTHETIC", "external_finemapping.pdf"), warn = FALSE, skipNul = TRUE))
vector_svg <- readLines(file.path(custom_out, "SYNTHETIC", "external_finemapping.svg"), warn = FALSE)
assert(any(grepl("/MediaBox [0 0 864 576]", vector_pdf, fixed = TRUE, useBytes = TRUE)) && any(grepl('width="864(pt)?" height="576(pt)?"', vector_svg)), "Vector dimensions must be pixels / DPI in inches")
custom_text <- readLines(file.path(custom_out, "SYNTHETIC", "INTERPRETATION.txt"))
assert(any(grepl("80% credible sets", custom_text, fixed = TRUE)) && !any(grepl("95% credible sets", custom_text, fixed = TRUE)), "Credible-set interpretation retained stale coverage")
assert(any(grepl("beta +/- 1.959964 * SE", custom_text, fixed = TRUE)), "Changing CS coverage altered external beta interval semantics")
assert(all(vapply(c("PNG", "PDF", "SVG"), function(x) any(grepl(paste0("[", x, "]("), readLines(file.path(custom_out, "index.md")), fixed = TRUE)), logical(1))), "Index omitted selected formats")
# CWD-relative overrides are distinct from config-relative defaults, including the source() helper.
cwd <- getwd(); setwd(folder)
helper <- run_external_finemapping_config(config_path, overrides = list(out = "helper_relative", plot_formats = "pdf", max_effects = 1L))
setwd(cwd)
assert(helper$effective_settings$paths$out == normalizePath(file.path(folder, "helper_relative")), "R helper override path was incorrectly relative to config directory")
assert(file.exists(file.path(folder, "helper_relative", "SYNTHETIC", "external_finemapping.pdf")) && !file.exists(file.path(folder, "helper_relative", "SYNTHETIC", "external_finemapping.png")), "Unrequested PNG was manufactured")
assert(!any(grepl("[PNG]", readLines(file.path(folder, "helper_relative", "index.md")), fixed = TRUE)), "Index links to unrequested PNG")

# Selection uses the matching source rows/columns, including a freshly correct eigendecomposition.
good <- load_fixture(write_fixture())
for (selection in c("snps", "indels")) {
  selected <- external_select_variants(good, selection)
  indices <- which(if (selection == "snps") good$variants$variant_type == "SNP" else good$variants$variant_type != "SNP")
  assert(identical(selected$variants$variant_id, good$variants$variant_id[indices]) && max(abs(selected$R-R[indices, indices])) < 1e-12, "Variant class subset broke LD alignment")
  decomposition <- attr(selected$R, "eigen")
  assert(max(abs(decomposition$vectors %*% diag(decomposition$values) %*% t(decomposition$vectors)-selected$R)) < 1e-10, "Subset reused stale full-LD eigenvectors")
  subset_out <- file.path(folder, selection)
  result <- run_external_finemapping(input, subset_out, 1, variant_types = selection, max_effects = 2)
  subset_pip <- read.delim(file.path(subset_out, "SYNTHETIC", "variant_pip.tsv"))
  subset_summary <- jsonlite::fromJSON(file.path(subset_out, "SYNTHETIC", "summary.json"))
  assert(identical(subset_pip$variant_id, good$variants$variant_id[indices]) && subset_summary$n_source_variants == 8 && subset_summary$n_variants == length(indices), "Source/selected counts or class rows incorrect")
  assert(identical(as.integer(subset_summary$selected_source_variant_indices), indices), "Source index provenance missing")
  subset_qc <- EXTERNAL_QC; subset_qc$max_effects <- 2L
  expected <- external_usable_fit(external_fit(selected$variants$beta/selected$variants$se, selected$R, m$effective_n, subset_qc), length(indices))
  assert(all(is.finite(subset_pip$pip)) && max(abs(expected-subset_pip$pip)) < 1e-12, "Class-specific fit used incorrect matrix or variant ordering")
}
# Neither an invalid omitted variant nor invalid cross-class LD can escape full validation.
bad <- variants; bad$se[2] <- 0; write_fixture(bad)
subset_bad <- file.path(folder, "subset_bad_variant"); run_external_finemapping(input, subset_bad, 1, variant_types = "snps")
assert(jsonlite::fromJSON(file.path(subset_bad, "SYNTHETIC", "summary.json"))$status == "unavailable_input_validation", "Class selection bypassed invalid omitted variant")
bad_R <- R; bad_R[1, 2] <- NA; write_fixture(r = bad_R)
subset_bad_ld <- file.path(folder, "subset_bad_ld"); run_external_finemapping(input, subset_bad_ld, 1, variant_types = "snps")
assert(jsonlite::fromJSON(file.path(subset_bad_ld, "SYNTHETIC", "summary.json"))$status == "unavailable_input_validation", "Class selection bypassed invalid cross-class LD")
# Zero or one selected variant remains explicitly unavailable, with any available external beta preserved.
for (n_selected in c(0L, 1L)) {
  few <- variants; few$ref[2:3] <- "A"; few$alt[2:3] <- "T"; few$variant_type[2:3] <- "SNP"
  if (n_selected == 1L) { few$alt[2] <- "AT"; few$variant_type[2] <- "INDEL" }
  write_fixture(few); few_out <- file.path(folder, paste0("few_", n_selected))
  run_external_finemapping(input, few_out, 1, variant_types = "indels", plot_formats = c("png", "pdf"))
  few_summary <- jsonlite::fromJSON(file.path(few_out, "SYNTHETIC", "summary.json"))
  few_pip <- read.delim(file.path(few_out, "SYNTHETIC", "variant_pip.tsv"))
  assert(few_summary$status == "unavailable_too_few_selected_variants" && few_summary$n_source_variants == 8 && few_summary$n_variants == n_selected, "Insufficient selected variants mislabeled")
  assert(nrow(few_pip) == n_selected && all(is.na(few_pip$pip)) && !file.exists(file.path(few_out, "SYNTHETIC", "model_effective_n.rds")), "Insufficient selected variants manufactured PIPs")
  if (n_selected) assert(few_pip$beta == variants$beta[2], "Unavailable model erased an observed external beta")
  assert(file.info(file.path(few_out, "SYNTHETIC", "external_finemapping.png"))$size > 1000 && grepl("Insufficient selected variants", few_summary$reason, fixed = TRUE), "Insufficient-selection figure/reason missing")
}
write_fixture()
# max_iter reaches the actual call: forced one-iteration nonconvergence cannot expose internal PIPs.
iteration_out <- file.path(folder, "iteration_limit")
run_external_finemapping(input, iteration_out, 1, max_iter = 1)
iteration_summary <- jsonlite::fromJSON(file.path(iteration_out, "SYNTHETIC", "summary.json"))
assert(iteration_summary$status == "unavailable_fit" && iteration_summary$max_iter == 1 && all(is.na(read.delim(file.path(iteration_out, "SYNTHETIC", "variant_pip.tsv"))$pip)), "Requested iteration limit not applied to convergence gate")

# Strict option/schema validation rejects typos, unsafe types and attempts to configure fixed gates.
base <- list(input_dir = input, out = file.path(folder, "unused"))
invalid_options <- list(list(workers = 0), list(workers = TRUE), list(workers = 1.5), list(max_effects = 0), list(max_effects = Inf),
  list(coverage = 0), list(coverage = 1), list(coverage = NA_real_), list(coverage = "0.9"), list(max_iter = -1),
  list(plot_width = 0), list(plot_height = c(600, 700)), list(plot_dpi = Inf), list(plot_formats = c("png", "png")),
  list(plot_formats = "jpg"), list(plot_formats = character()), list(variant_types = "snp"), list(locus = "../unsafe"), list(s_block = 1))
for (bad_options in invalid_options) fails(external_resolve_options(c(base, bad_options)), "Invalid setting accepted")
bad_configs <- list(list(qc = list(s_block = 1)), list(analysis = list(residual_variance = TRUE)), list(plot = list(formats = "png")),
  list(plot = list(formats = list(named = "png"))), list(paths = list(input_dir = NULL)), list(analysis = list(workers = "2")))
bad_config_path <- file.path(config_dir, "invalid.json")
for (bad_config in bad_configs) {
  external_json(bad_config, bad_config_path)
  fails(external_resolve_options(base, bad_config_path), "Invalid JSON setting accepted")
}
writeLines('{"analysis":{"coverage":0.8,"coverage":0.9}}', bad_config_path)
fails(external_resolve_options(base, bad_config_path), "Duplicate JSON keys accepted")
fails(external_cli(c("--out", "unused", "--disable-qc", "true")), "Unknown CLI QC bypass accepted")
fails(external_cli(c("--workers", "1", "--workers", "2")), "Duplicate CLI argument accepted")
fails(external_cli(c("--plot-formats", "png,")), "Empty trailing format accepted")
fails(external_cli("--input-dir"), "Missing argument value accepted")
fails(run_external_finemapping_config(config_path, overrides = list(qc = list(s_block = 1))), "R helper QC bypass accepted")
# Configuration mutation leaves an incomplete run, with no misleading completion marker.
external_fit <- function(...) { result <- real_fit(...); cat(" ", file = config_path, append = TRUE); result }
mutation_out <- file.path(folder, "configuration_changed")
fails(run_external_finemapping_config(config_path, list(out = mutation_out)), "Changed configuration was certified complete", "Configuration changed during analysis")
external_fit <- real_fit
assert(!file.exists(file.path(mutation_out, "COMPLETED.txt")), "Changed configuration retained completion marker")
# Help is standalone and sourcing above did not execute the CLI.
help_output <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(EXTERNAL_SCRIPT_PATH), "--help"), stdout = TRUE, stderr = TRUE)
assert(is.null(attr(help_output, "status")) && all(vapply(c("--config", "--max-effects", "--variant-types", "--plot-formats", "RStudio"), function(x) any(grepl(x, help_output, fixed = TRUE)), logical(1))), "Standalone CLI help failed or omitted key options")
on.exit_cleanup()
cat("PASS: joint SNP/indel SuSiE; fixed provenance/LD gates; configuration/default parity; CLI/RStudio overrides; exact config snapshots; selected LD alignment; insufficient variants; configurable models and PNG/PDF/SVG; locally deterministic diagnostics with caller RNG restoration; preserved runs.\n")
