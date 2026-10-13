#!/usr/bin/env Rscript
# Small native R fixtures for the evidence boundary and public provenance.
entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(entry)), "audit_finemapping.R"))

assert <- function(value, message) if (!isTRUE(value)) stop(message)
expect_error <- function(expression) {
  failed <- tryCatch({ force(expression); FALSE }, error = function(e) TRUE)
  assert(failed, "Expected invalid input to be rejected")
}

temporary <- tempfile("audit_fixture_"); dir.create(temporary)
on.exit_cleanup <- function() unlink(temporary, recursive = TRUE)
run <- file.path(temporary, "cohort"); dir.create(run)
legacy <- file.path(temporary, "legacy"); dir.create(legacy)
locus <- file.path(legacy, "loci", "chr1_12"); dir.create(locus, recursive = TRUE)
out <- file.path(temporary, "audit")
gwas <- file.path(temporary, "gwas.tsv")

summary <- list(synthetic = FALSE, inputs_stable_during_read = TRUE,
  genome_build_user_declared = "GRCh38", group_counts = list(T2D = 2L, Prediabetes = 2L, Control = 2L),
  haplotype_carrier_tests_performed = 0L, scan_counts = list(high_quality_homozygous_reference_genotypes = 0L))
save_summary <- function(x = summary) jsonlite::write_json(x, file.path(run, "run_summary.json"), auto_unbox = TRUE)
save_summary(); writeLines("Completed", file.path(run, "COMPLETED.txt"))
lead <- list(locus_id = "L0001", chrom = "chr1", lead_pos = "12", lead_effect_allele = "A",
  lead_other_allele = "G", lead_p = "3.35E-1003", status = "lead_below_min_call_rate", selected_markers = "0")
write_rows(file.path(run, "lead_loci.tsv"), list(lead), names(lead))
gwas_row <- list(Chromsome = "1", Position = "12", EffectAllele = "A", NonEffectAllele = "G",
  Pval = "3.35e-1003", Beta = "-0.2", SE = "0.01", Neff = "12345", liftover_strand = "+")
write_rows(gwas, list(gwas_row), names(gwas_row))
legacy_row <- list(CHR = "chr1", LEAD_BP = "12", STATUS = "SUCCESS")
write_rows(file.path(legacy, "all_loci_summary.tsv"), list(legacy_row), names(legacy_row))
qc <- list(ONT_samples = "6", Variants_analyzed = "1", Converged = "TRUE", GWAS_n = "NOT_SUPPLIED", LD_shrinkage = "0.02")
write_rows(file.path(locus, "susie_QC.tsv"), list(qc), names(qc))
pip <- list(Variant_key = "chr1:12:A:G", GWAS_EA = "A", GWAS_NEA = "G", ONT_REF = "A", ONT_ALT = "G",
  Beta = "-0.2", SE = "0.01", Z_ALT = "20", PIP = "0.99999123456789")
write_rows(file.path(locus, "susie_PIP.tsv"), list(pip), names(pip))
geno <- file.path(locus, "locus_genotypes.tsv")
writeLines("chr1\t12\tA\tG\t0|1\t1/0\t1|1\t./.\t.|.\t1/1", geno)
genotype_hash <- file_sha256(geno)
old_hap <- file.path(temporary, "Haplotype_comparison"); dir.create(old_hap)
writeLines(c('risk_code_anchor <- if (g$Effect_allele_in_ONT[anchor_i] == "ALT") 1 else 0',
  "pval_num <- as.numeric(g$Pval)"), file.path(old_hap, "run_one_locus_haplotype.R"))

# Known digest detects accidental MD5 substitution or hash/string conversion.
abc <- file.path(temporary, "abc"); writeBin(charToRaw("abc"), abc)
assert(identical(file_sha256(abc), "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"), "SHA256 digest incorrect")

gt <- legacy_genotypes(geno, "chr1:12:A:G", 6L)
assert(gt$ref_hom == 0L && gt$het == 2L && gt$alt_hom == 2L && gt$missing_or_other == 2L,
  "Missing genotypes were misclassified as reference")
assert(legacy_genotypes(geno, "chr1:12:A:G", 7L)$width_mismatch_rows == 1L, "Sample-count mismatch was missed")
assert(legacy_alignment(read_tsv(file.path(locus, "susie_PIP.tsv")))$errors == 0L, "Correct REF-oriented z failed")
wrong <- read_tsv(file.path(locus, "susie_PIP.tsv")); wrong$Z_ALT <- "-20"
assert(legacy_alignment(wrong)$errors == 1L, "Allele flip error was missed")

result <- audit_finemapping(run, gwas, legacy, out)
assert(result$resolved_candidate_variants == 0L && result$legacy_eligible_for_causal_ranking == 0L,
  "High unvalidated PIP or convergence became a causal candidate")
assert(nrow(read_tsv(file.path(out, "candidate_variants.tsv"))) == 0L, "Candidate table must be header-only")
assert(result$legacy_exclusion_counts$no_reference_homozygotes_in_retained_ld_genotypes == 1L, "Reference-call absence lost")
external <- read_tsv(file.path(out, "external_gwas_leads.tsv"))
assert(external$external_gwas_p == "3.35e-1003" && external$effect_direction == "negative", "Tiny P or Beta sign lost")
assert(identical(file_sha256(geno), genotype_hash), "Audit modified original genotypes")
assert(identical(names(result$provenance), c("cohort_summary", "selected_leads", "gwas")), "Provenance role schema differs")
assert(is.null(result$provenance_md5), "Absolute-path MD5 provenance was retained")
for (name in names(result$output_sha256)) assert(identical(result$output_sha256[[name]], file_sha256(file.path(out, name))), "Output checksum mismatch")
assert(isTRUE(result$legacy_haplotype_script_review$effect_allele_assigned_risk_without_beta_sign), "Historical risk-label caveat lost")
assert(length(result$method_sources) == 2L, "Primary method links missing")
for (path in list.files(out, full.names = TRUE)) {
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  assert(!grepl(temporary, text, fixed = TRUE), "Private filesystem path exported")
  assert(!grepl(pip$PIP, text, fixed = TRUE), "Unvalidated numeric PIP exported")
}
expect_error(audit_finemapping(run, gwas, legacy, out))

# Invalid sources fail before any new output folder is created.
fresh <- file.path(temporary, "fresh")
bad <- summary; bad$synthetic <- TRUE; save_summary(bad)
expect_error(audit_finemapping(run, gwas, legacy, fresh)); assert(!dir.exists(fresh), "Synthetic output was created")
bad <- summary; bad$haplotype_carrier_tests_performed <- 1L; save_summary(bad)
expect_error(audit_finemapping(run, gwas, legacy, fresh)); assert(!dir.exists(fresh), "Tested cohort was mislabeled unresolved")
save_summary()
changed <- gwas_row; changed$Pval <- "1e-8"; write_rows(gwas, list(changed), names(changed))
expect_error(audit_finemapping(run, gwas, legacy, fresh)); assert(!dir.exists(fresh), "Mismatched GWAS output was created")
write_rows(gwas, list(gwas_row, gwas_row), names(gwas_row))
expect_error(audit_finemapping(run, gwas, legacy, fresh)); assert(!dir.exists(fresh), "Duplicate GWAS output was created")
on.exit_cleanup()
cat("PASS: native R fine-mapping audit, tiny-P/allele checks, missing-call handling, causal boundaries, SHA256 and public provenance.\n")
