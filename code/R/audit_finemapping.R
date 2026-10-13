#!/usr/bin/env Rscript
# Independent R audit of saved GWAS and legacy fine-mapping artifacts.
# Never recomputes or promotes unvalidated legacy PIPs to causal evidence.
entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
if (!exists("read_tsv", mode = "function")) source(file.path(dirname(normalizePath(entry)), "compare_haplotypes.R"))

LEGACY_FIELDS <- strsplit("legacy_locus run_status qc_present model_converged ld_sample_size gwas_n ld_shrinkage variants_analyzed pip_rows_present credible_set_memberships_present retained_genotype_rows_found retained_gt_ref_hom retained_gt_het retained_gt_alt_hom retained_gt_missing_or_other retained_gt_width_mismatch_rows allele_or_z_alignment_errors palindromic_retained_variants eligible_for_causal_ranking exclusion_reasons", " ")[[1]]
LEAD_FIELDS <- strsplit("locus_id chrom position effect_allele other_allele external_gwas_p external_gwas_beta external_gwas_se external_gwas_eaf external_gwas_ncases external_gwas_ncontrols external_gwas_neff effect_direction liftover_strand allele_orientation_review cohort_locus_status cohort_association_evidence causal_status", " ")[[1]]
CANDIDATE_FIELDS <- c("locus_id", "variant", "evidence", "validated_pip", "causal_status")

METHOD_SOURCES <- list(
  list(title = "SuSiE-RSS summary-statistic and LD diagnostics",
    url = "https://stephenslab.github.io/susieR/articles/susierss_diagnostic.html",
    application = "Validate GWAS allele coding and reference LD consistency before interpreting fine-mapping."),
  list(title = "SuSiE-RSS reference documentation",
    url = "https://stephenslab.github.io/susieR/reference/susie_rss.html",
    application = "GWAS sample size is recommended, although optional; omission alone does not invalidate a model."))

file_sha256 <- function(path) {
  if (requireNamespace("openssl", quietly = TRUE)) {
    con <- file(path, "rb"); on.exit(close(con))
    return(unname(unclass(as.character(openssl::sha256(con)))))
  }
  command <- Sys.which("sha256sum")
  if (!nzchar(command)) stop("Audit provenance requires R package openssl or the sha256sum command.")
  output <- system2(command, shQuote(path), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(output, "status")) || !length(output) || !grepl("^[0-9a-f]{64} ", output[1]))
    stop("Unable to calculate SHA256 provenance.")
  substr(output[1], 1L, 64L)
}

file_fingerprint <- function(path) list(name = basename(path), bytes = unname(file.info(path)$size), sha256 = file_sha256(path))

legacy_method_metadata <- function(folder) {
  summary <- file.path(folder, "all_loci_summary.tsv")
  if (!file.exists(summary)) return(list(available = FALSE, reason = "legacy_summary_not_supplied"))
  script <- file.path(folder, "run_one_locus.R")
  text <- if (file.exists(script)) paste(readLines(script, warn = FALSE), collapse = "\n") else ""
  list(available = TRUE, summary_source = file_fingerprint(summary),
    script_source = if (file.exists(script)) file_fingerprint(script) else NULL,
    script_snapshot_is_not_proof_of_historical_execution = TRUE,
    script_mean_imputes_missing_dosages = grepl("dosage[i,missing] <- mean(", text, fixed = TRUE),
    script_has_dp_gq_filters = grepl("min_dp|min_gq", text),
    script_explicitly_checks_ld_z_consistency = grepl("estimate_s_rss\\(|kriging_rss\\(", text),
    script_mentions_liftover_strand_or_palindromic_checks = grepl("liftover_strand|palindrom", text))
}

historical_haplotype_metadata <- function(folder) {
  script <- file.path(dirname(folder), "Haplotype_comparison", "run_one_locus_haplotype.R")
  if (!file.exists(script)) return(list(available = FALSE))
  text <- paste(readLines(script, warn = FALSE), collapse = "\n")
  list(available = TRUE, source = file_fingerprint(script),
    effect_allele_assigned_risk_without_beta_sign =
      grepl('risk_code_anchor <- if (g$Effect_allele_in_ONT[anchor_i] == "ALT") 1 else 0', text, fixed = TRUE) &&
      !grepl("g$Beta", text, fixed = TRUE),
    gwas_p_converted_to_machine_numeric = grepl("as.numeric(g$Pval)", text, fixed = TRUE),
    phase_set_continuity_validation_present = grepl("phase_set", text, fixed = TRUE), results_used_in_current_analysis = FALSE)
}

complete_audit_metadata <- function(result, run_dir, gwas, legacy_dir, out) {
  result$schema_version <- 1L
  result$provenance_md5 <- NULL
  result$provenance <- list(cohort_summary = file_fingerprint(file.path(run_dir, "run_summary.json")),
    selected_leads = file_fingerprint(file.path(run_dir, "lead_loci.tsv")), gwas = file_fingerprint(gwas))
  result$output_sha256 <- setNames(lapply(c("external_gwas_leads.tsv", "legacy_locus_audit.tsv", "candidate_variants.tsv"),
    function(name) file_sha256(file.path(out, name))), c("external_gwas_leads.tsv", "legacy_locus_audit.tsv", "candidate_variants.tsv"))
  result$legacy_method <- legacy_method_metadata(legacy_dir)
  result$legacy_haplotype_script_review <- historical_haplotype_metadata(legacy_dir)
  result$method_sources <- METHOD_SOURCES
  if (isTRUE(result$legacy_haplotype_script_review$effect_allele_assigned_risk_without_beta_sign))
    result$limitations <- unique(c(result$limitations,
      "Historical haplotype code labels the effect allele as risk without Beta sign, converts tiny p-values to machine numeric, and lacks phase-set continuity validation. Its group comparisons are not reused as current evidence."))
  result
}

write_audit_report <- function(result, out) {
  jsonlite::write_json(result, file.path(out, "audit_summary.json"), pretty = TRUE, auto_unbox = TRUE, na = "null", null = "null")
  writeLines(c("Independent R variant-evidence audit", "", "No causal variants can be resolved from the validated available evidence.",
    sprintf("%d external GWAS leads; %d legacy loci audited; 0 eligible for causal ranking.", result$external_gwas_leads_reported, result$legacy_loci),
    "candidate_variants.tsv is intentionally header-only; this means unresolved, not no causal variants.", "", result$limitations,
    "", "Method sources:", vapply(METHOD_SOURCES, function(s) paste(s$title, s$url, sep = ": "), character(1))), file.path(out, "REPORT.txt"))
}

normalize_chrom <- function(x) paste0("chr", sub("^chr", "", x))
is_palindromic <- function(a, b) paste(sort(c(a, b)), collapse = "") %in% c("AT", "CG")

external_lead_evidence <- function(leads, gwas) {
  x <- read_tsv(gwas)
  if (!all(c("Chromsome", "Position", "EffectAllele", "NonEffectAllele", "Pval", "Beta", "SE") %in% names(x))) stop("Audit requires original GWAS column names.")
  key <- function(chrom, pos, a, b) paste(normalize_chrom(chrom), pos, toupper(a), toupper(b), sep = ":")
  source_key <- key(x$Chromsome, x$Position, x$EffectAllele, x$NonEffectAllele)
  wanted <- key(leads$chrom, leads$lead_pos, leads$lead_effect_allele, leads$lead_other_allele)
  if (anyDuplicated(source_key[source_key %in% wanted])) stop("Ambiguous duplicate selected GWAS lead.")
  ix <- match(wanted, source_key)
  if (anyNA(ix)) stop("Selected lead lacks an exact GWAS allele match.")
  lapply(seq_len(nrow(leads)), function(i) {
    l <- as.list(leads[i, ]); s <- as.list(x[ix[i], ])
    if (!identical(decimal_parts(s$Pval), decimal_parts(l$lead_p))) stop("Lead and GWAS p-values differ.")
    beta <- as.numeric(s$Beta); se <- as.numeric(s$SE)
    if (!is.finite(beta) || !is.finite(se) || se <= 0) stop("Invalid GWAS effect or uncertainty.")
    strand <- s$liftover_strand %||% "NA"
    orientation <- if (is_palindromic(l$lead_effect_allele, l$lead_other_allele)) "palindromic_requires_forward_reference_validation"
      else if (strand == "-") "negative_liftover_strand_requires_harmonization" else "exact_allele_match_only_not_reference_attestation"
    list(locus_id = l$locus_id, chrom = normalize_chrom(l$chrom), position = l$lead_pos, effect_allele = l$lead_effect_allele,
      other_allele = l$lead_other_allele, external_gwas_p = s$Pval, external_gwas_beta = s$Beta, external_gwas_se = s$SE,
      external_gwas_eaf = s$EAF %||% "NA", external_gwas_ncases = s$Ncases %||% "NA", external_gwas_ncontrols = s$Ncontrols %||% "NA",
      external_gwas_neff = s$Neff %||% "NA", effect_direction = if (beta > 0) "positive" else if (beta < 0) "negative" else "zero",
      liftover_strand = strand, allele_orientation_review = orientation, cohort_locus_status = l$status,
      cohort_association_evidence = "not_estimated", causal_status = "unresolved_external_association_only")
  })
}

legacy_alignment <- function(pip) {
  errors <- palindromic <- 0L
  if (!nrow(pip)) return(list(errors = errors, palindromic = palindromic))
  for (i in seq_len(nrow(pip))) {
    r <- as.list(pip[i, ])
    ea <- toupper(r$GWAS_EA %||% ""); nea <- toupper(r$GWAS_NEA %||% "")
    ref <- toupper(r$ONT_REF %||% ""); alt <- toupper(r$ONT_ALT %||% "")
    if (is_palindromic(ref, alt)) palindromic <- palindromic+1L
    beta <- suppressWarnings(as.numeric(r$Beta %||% NA_character_)); se <- suppressWarnings(as.numeric(r$SE %||% NA_character_))
    z <- suppressWarnings(as.numeric(r$Z_ALT %||% NA_character_)); expected <- beta/se * if (ea == alt) 1 else -1
    ok <- identical(sort(c(ea, nea)), sort(c(ref, alt))) && ref != alt && is.finite(se) && se > 0 &&
      is.finite(z) && is.finite(expected) && abs(z-expected) <= max(1e-8, 1e-6*max(abs(z), abs(expected)))
    if (!ok) errors <- errors+1L
  }
  list(errors = errors, palindromic = palindromic)
}

legacy_genotypes <- function(path, keys, expected_n) {
  if (!file.exists(path)) return(list(available = FALSE))
  values <- setNames(rep(0L, 6L), c("rows", "ref_hom", "het", "alt_hom", "missing_or_other", "width_mismatch_rows"))
  con <- file(path, "rt"); on.exit(close(con))
  seen <- character()
  repeat {
    lines <- readLines(con, n = 1000L, warn = FALSE)
    if (!length(lines)) break
    for (line in lines) {
      f <- strsplit(line, "\t", fixed = TRUE)[[1]]
      if (length(f) < 5L) stop("Malformed legacy genotype row.")
      alleles <- sort(toupper(f[3:4])); key <- paste(f[1], f[2], alleles[1], alleles[2], sep = ":")
      if (!key %in% keys) next
      if (key %in% seen) stop("Duplicate retained legacy genotype variant.")
      seen <- c(seen, key); values["rows"] <- values["rows"]+1L
      if (!is.na(expected_n) && length(f)-4L != expected_n) values["width_mismatch_rows"] <- values["width_mismatch_rows"]+1L
      gt <- gsub("|", "/", f[-(1:4)], fixed = TRUE)
      types <- c(`0/0` = "ref_hom", `0/1` = "het", `1/0` = "het", `1/1` = "alt_hom")[gt]
      types[is.na(types)] <- "missing_or_other"
      counts <- table(types)
      for (k in names(counts)) values[k] <- values[k]+counts[[k]]
    }
  }
  c(list(available = TRUE), as.list(values))
}

audit_legacy <- function(folder) {
  path <- file.path(folder, "all_loci_summary.tsv")
  if (!file.exists(path)) return(list(rows = list(), method = list(available = FALSE, reason = "legacy_summary_not_supplied")))
  summary <- read_tsv(path); rows <- list()
  method <- legacy_method_metadata(folder)
  for (i in seq_len(nrow(summary))) {
    record <- as.list(summary[i, ]); name <- paste(record$CHR, record$LEAD_BP, sep = "_")
    locus <- file.path(folder, "loci", name)
    read_optional <- function(name) if (file.exists(file.path(locus, name))) read_tsv(file.path(locus, name)) else data.frame()
    qc_rows <- read_optional("susie_QC.tsv")
    if (nrow(qc_rows) > 1L) stop("Expected one QC row per locus.")
    qc <- if (nrow(qc_rows)) as.list(qc_rows[1, ]) else list()
    pip <- read_optional("susie_PIP.tsv"); cs <- read_optional("susie_credible_sets.tsv")
    keys <- pip$Variant_key %||% character()
    if (anyDuplicated(keys)) stop("Duplicate PIP variant keys.")
    alignment <- legacy_alignment(pip)
    n <- suppressWarnings(as.integer(qc$ONT_samples %||% NA_character_))
    gt <- legacy_genotypes(file.path(locus, "locus_genotypes.tsv"), keys, n)
    reasons <- c("ld_callability_and_gwas_population_match_unvalidated", "ld_z_consistency_not_validated", "functional_causal_evidence_not_supplied")
    if (record$STATUS != "SUCCESS") reasons <- c(reasons, "legacy_fit_failed")
    if (!length(qc)) reasons <- c(reasons, "legacy_qc_missing")
    else if ((qc$Converged %||% "") != "TRUE") reasons <- c(reasons, "model_not_converged")
    if (identical(qc$GWAS_n, "NOT_SUPPLIED")) reasons <- c(reasons, "gwas_sample_size_not_supplied")
    if (length(keys) && !gt$available) reasons <- c(reasons, "ld_genotype_artifact_missing")
    else if (length(keys)) {
      if (gt$rows != length(keys)) reasons <- c(reasons, "retained_ld_variants_missing_from_genotypes")
      if (gt$width_mismatch_rows > 0) reasons <- c(reasons, "ld_sample_count_mismatch")
      if (gt$ref_hom == 0) reasons <- c(reasons, "no_reference_homozygotes_in_retained_ld_genotypes")
    }
    if (alignment$errors > 0) reasons <- c(reasons, "allele_or_z_alignment_error")
    if (alignment$palindromic > 0) reasons <- c(reasons, "palindromic_orientation_unvalidated")
    rows[[i]] <- list(legacy_locus = name, run_status = record$STATUS, qc_present = length(qc) > 0,
      model_converged = qc$Converged %||% NA_character_, ld_sample_size = n, gwas_n = qc$GWAS_n %||% NA_character_,
      ld_shrinkage = qc$LD_shrinkage %||% NA_character_, variants_analyzed = qc$Variants_analyzed %||% NA_character_,
      pip_rows_present = nrow(pip), credible_set_memberships_present = nrow(cs), retained_genotype_rows_found = gt$rows %||% NA_integer_,
      retained_gt_ref_hom = gt$ref_hom %||% NA_integer_, retained_gt_het = gt$het %||% NA_integer_, retained_gt_alt_hom = gt$alt_hom %||% NA_integer_,
      retained_gt_missing_or_other = gt$missing_or_other %||% NA_integer_, retained_gt_width_mismatch_rows = gt$width_mismatch_rows %||% NA_integer_,
      allele_or_z_alignment_errors = alignment$errors, palindromic_retained_variants = alignment$palindromic,
      eligible_for_causal_ranking = FALSE, exclusion_reasons = paste(reasons, collapse = ";"))
  }
  list(rows = rows, method = method)
}

audit_finemapping <- function(run_dir, gwas, legacy_dir, out) {
  if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) stop("Use a new or empty audit output directory.")
  if (!file.exists(file.path(run_dir, "COMPLETED.txt"))) stop("Incomplete cohort run.")
  source <- jsonlite::fromJSON(file.path(run_dir, "run_summary.json"), simplifyVector = FALSE)
  leads <- read_tsv(file.path(run_dir, "lead_loci.tsv"))
  if (isTRUE(source$synthetic) || grepl("demo", basename(run_dir), ignore.case = TRUE)) stop("Synthetic results cannot supply cohort causal evidence.")
  if ((source$haplotype_carrier_tests_performed %||% 0) > 0 || any(as.numeric(leads$selected_markers) > 0))
    stop("This unresolved-cohort audit cannot evaluate tested haplotypes; independent evidence evaluation is required.")
  external <- external_lead_evidence(leads, gwas)
  legacy <- audit_legacy(legacy_dir)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write_rows(file.path(out, "external_gwas_leads.tsv"), external, LEAD_FIELDS)
  write_rows(file.path(out, "legacy_locus_audit.tsv"), legacy$rows, LEGACY_FIELDS)
  write_rows(file.path(out, "candidate_variants.tsv"), list(), CANDIDATE_FIELDS)
  result <- list(implementation = "independent_R", created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    analysis = "aggregate_evidence_audit_not_new_finemapping", synthetic = FALSE, cohort_group_counts = source$group_counts,
    cohort_selected_loci = nrow(leads), cohort_haplotype_association_tests = 0, cohort_locus_status_counts = counts_list(leads$status),
    cohort_high_quality_reference_homozygotes = source$scan_counts$high_quality_homozygous_reference_genotypes,
    external_gwas_leads_reported = length(external), resolved_candidate_variants = 0,
    candidate_status = "not_resolvable_from_validated_available_evidence", legacy_loci = length(legacy$rows),
    legacy_status_counts = counts_list(vapply(legacy$rows, `[[`, character(1), "run_status")), legacy_eligible_for_causal_ranking = 0,
    legacy_exclusion_counts = counts_list(unlist(strsplit(vapply(legacy$rows, `[[`, character(1), "exclusion_reasons"), ";", fixed = TRUE))),
    legacy_method = legacy$method,
    limitations = c("No eligible cohort haplotypes: no group difference or causal variant can be resolved.",
      "GWAS leads are external associations; phenotype coding/effect scale and reference orientation require validation.",
      "Legacy convergence/PIPs cannot resolve unvalidated callability, population matching, or LD/z consistency.",
      "No reference homozygotes is a callability concern; missing genotypes must not be replaced with reference calls.",
      "Missing GWAS sample size is a limitation, not by itself proof of invalid fine-mapping.",
      "Legacy variants/genotypes are locus memberships: overlapping loci may count a variant more than once."))
  result$genome_build_declared_in_cohort_run <- source$genome_build_user_declared %||% "not_supplied"
  result <- complete_audit_metadata(result, run_dir, gwas, legacy_dir, out)
  write_audit_report(result, out)
  cat("Independent R audit: ", length(legacy$rows), " legacy loci; 0 eligible causal candidates.\n", sep = "")
  invisible(result)
}

if (sys.nframe() == 0L) {
  argv <- commandArgs(TRUE)
  if ("--help" %in% argv) {
    cat("Rscript code/R/audit_finemapping.R --run-dir RUN --gwas GWAS_TSV --legacy-dir LEGACY_FOLDER --out NEW_DIRECTORY\n")
  } else {
    if (length(argv) != 8L || !setequal(argv[seq(1, 8, 2)], c("--run-dir", "--gwas", "--legacy-dir", "--out"))) stop("Use --help for required arguments.")
    args <- as.list(argv[seq(2, 8, 2)]); names(args) <- gsub("-", "_", sub("^--", "", argv[seq(1, 8, 2)]))
    do.call(audit_finemapping, args)
  }
}
