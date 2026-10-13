#!/usr/bin/env Rscript
# Exploratory, window-limited external GWAS SuSiE-RSS with independent ALT-coded LD.
# No participant-carrier associations or cohort genotype files are modified.
EXTERNAL_SCRIPT_PATH <- if (sys.nframe() > 0L && !is.null(sys.frame(1)$ofile)) {
  normalizePath(sys.frame(1)$ofile, mustWork = TRUE)
} else {
  normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = TRUE)
}
EXTERNAL_FIELDS <- c("variant_id", "chrom", "pos", "ref", "alt", "variant_type", "beta", "se", "pval", "af", "gwas_n", "reference_af")
EXTERNAL_QC <- list(min_maf = .01, max_af_difference = .1, matrix_tolerance = 1e-8,
  s_warning = .1, s_block = .5, kriging_abs_residual_warning = 6, kriging_logLR_warning = 2,
  kriging_abs_z_for_logLR = 2, boundary_distance_bp = 5000,
  strong_association_p = 5e-8, max_iter = 200L, coverage = .95)
EXTERNAL_METHODS <- c("https://stephenslab.github.io/susieR/reference/susie_rss.html",
  "https://stephenslab.github.io/susieR/articles/susierss_diagnostic.html",
  "https://elifesciences.org/articles/88768")

external_sha256 <- function(path) {
  if (!requireNamespace("openssl", quietly = TRUE)) stop("R package openssl is required for SHA256 provenance.")
  con <- file(path, "rb"); on.exit(close(con))
  unname(unclass(as.character(openssl::sha256(con))))
}
external_write_table <- function(data, path) write.table(data, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
external_json <- function(data, path) jsonlite::write_json(data, path, pretty = TRUE, auto_unbox = TRUE, na = "null", null = "null", digits = 16)
external_capture <- function(expression) {
  warnings <- character()
  printed <- capture.output(value <- tryCatch(withCallingHandlers(force(expression), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }, message = function(m) {
    warnings <<- c(warnings, conditionMessage(m)); invokeRestart("muffleMessage")
  }), error = identity))
  list(value = value, warnings = unique(c(warnings, printed[nzchar(printed)])), error = if (inherits(value, "error")) conditionMessage(value) else NULL)
}
external_relative_file <- function(root, relative) {
  if (!is.character(relative) || length(relative) != 1L || is.na(relative) || startsWith(relative, "/")) stop("Input files must have relative paths.")
  path <- normalizePath(file.path(root, relative), mustWork = TRUE)
  if (!startsWith(path, paste0(normalizePath(root), "/"))) stop("Input file escapes input directory.")
  path
}
external_manifest <- function(input_dir) {
  path <- file.path(input_dir, "input_manifest.json")
  before <- external_sha256(path)
  manifest <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  if (external_sha256(path) != before) stop("Input manifest changed while reading.")
  expected <- list(genome_build = "GRCh38", effect_allele = "ALT", ld_dosage_allele = "ALT",
    allele_alignment = "exact_normalized_REF_ALT", reference_missingness = "complete_calls_only")
  for (field in names(expected)) if (!identical(manifest$provenance[[field]], expected[[field]])) stop("Missing or incompatible provenance: ", field)
  for (field in c("gwas_source", "ld_source")) {
    source <- manifest$provenance[[field]]
    valid_string <- is.character(source) && length(source) == 1L && !is.na(source) && nzchar(source)
    valid_record <- is.list(source) && length(source) > 0L && !is.null(names(source)) && all(nzchar(names(source)))
    if (!valid_string && !valid_record) stop("Missing source provenance: ", field)
  }
  counts <- unlist(manifest[c("n_cases", "n_controls")], use.names = FALSE)
  if (!is.numeric(counts) || length(counts) != 2L || any(!is.finite(counts) | counts <= 0 | counts != floor(counts))) stop("Positive integer external GWAS case and control counts are required.")
  if (!is.list(manifest$loci) || !length(manifest$loci)) stop("A nonempty loci manifest is required.")
  ids <- vapply(manifest$loci, function(x) x$locus_id, character(1))
  if (anyDuplicated(ids) || any(!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", ids)) || any(grepl("..", ids, fixed = TRUE))) stop("Locus identifiers must be unique and safe directory names.")
  list(data = manifest, sha256 = before, path = path,
       total_n = sum(counts), effective_n = 4/sum(1/counts))
}

external_load_locus <- function(input_dir, entry, total_n, qc = EXTERNAL_QC) {
  variant_file <- external_relative_file(input_dir, entry$variants_file)
  ld_file <- external_relative_file(input_dir, entry$ld_file)
  for (kind in c("variants", "ld")) {
    expected <- entry[[paste0(kind, "_sha256")]]
    path <- if (kind == "variants") variant_file else ld_file
    if (!is.character(expected) || length(expected) != 1L || !grepl("^[a-fA-F0-9]{64}$", expected) || external_sha256(path) != tolower(expected)) stop("SHA256 mismatch or missing hash: ", kind)
  }
  variants <- read.delim(variant_file, header = TRUE, check.names = FALSE, colClasses = "character", na.strings = c("NA", ".", ""))
  if (anyDuplicated(names(variants)) || !all(EXTERNAL_FIELDS %in% names(variants))) stop("Missing or duplicate variant table fields.")
  variants <- variants[, EXTERNAL_FIELDS, drop = FALSE]
  if (nrow(variants) < 2L) stop("At least two matched variants are required for locus fine-mapping and LD diagnostics.")
  if (anyNA(variants) || anyDuplicated(variants$variant_id) || any(!nzchar(variants$variant_id))) stop("Missing values or duplicate variant identifiers.")
  ld_ids <- unlist(entry$ld_variant_ids, use.names = FALSE)
  if (!identical(variants$variant_id, ld_ids)) stop("LD variant order must exactly match the variant table.")
  for (field in c("pos", "beta", "se", "pval", "af", "gwas_n", "reference_af")) {
    variants[[field]] <- suppressWarnings(as.numeric(variants[[field]]))
    if (any(!is.finite(variants[[field]]))) stop("Nonfinite or nonnumeric variant field: ", field)
  }
  if (length(unique(variants$chrom)) != 1L || any(variants$pos < 1 | variants$pos != floor(variants$pos))) stop("Invalid chromosome or position.")
  if (any(variants$se <= 0) || any(variants$pval < 0 | variants$pval > 1) || any(variants$gwas_n != total_n)) stop("Invalid SE, p-value or external GWAS sample size.")
  if (any(pmin(variants$af, 1-variants$af, variants$reference_af, 1-variants$reference_af) < qc$min_maf-1e-12)) stop("A variant fails the common-variant MAF requirement.")
  if (any(abs(variants$af-variants$reference_af) > qc$max_af_difference+1e-12)) stop("GWAS/reference ALT frequencies are incompatible.")
  if (any(!grepl("^[ACGT]+$", variants$ref)) || any(!grepl("^[ACGT]+$", variants$alt)) || any(variants$ref == variants$alt)) stop("Only explicit distinct REF/ALT DNA alleles are supported.")
  ref_len <- nchar(variants$ref); alt_len <- nchar(variants$alt)
  anchored <- substr(variants$ref, 1L, 1L) == substr(variants$alt, 1L, 1L) | substring(variants$ref, ref_len) == substring(variants$alt, alt_len)
  types <- ifelse(ref_len == 1L & alt_len == 1L, "SNP",
    ifelse(ref_len == 1L & alt_len > 1L & anchored, "insertion",
      ifelse(alt_len == 1L & ref_len > 1L & anchored, "deletion", ifelse(ref_len != alt_len, "complex_indel", "unsupported"))))
  supplied <- toupper(variants$variant_type)
  consistent <- (types == "SNP" & supplied == "SNP") | (types == "insertion" & supplied %in% c("INDEL", "INSERTION", "INS")) | (types == "deletion" & supplied %in% c("INDEL", "DELETION", "DEL")) | (types == "complex_indel" & supplied %in% c("INDEL", "COMPLEX_INDEL"))
  if (any(!consistent)) stop("Variant type disagrees with allele lengths or includes unsupported complex alleles.")
  if (anyDuplicated(paste(variants$chrom, variants$pos, variants$ref, variants$alt, sep = ":"))) stop("Duplicate normalized alleles with different identifiers.")
  variants$variant_type <- types
  window <- unlist(entry[c("window_start", "window_end")], use.names = FALSE)
  if (!is.numeric(window) || length(window) != 2L || any(!is.finite(window)) || window[1] < 1 || window[2] <= window[1] || any(window != floor(window)) || any(variants$pos < window[1] | variants$pos > window[2])) stop("Valid enclosing target-window coordinates are required.")
  con <- if (grepl("\\.gz$", ld_file)) gzfile(ld_file, "rt") else file(ld_file, "rt")
  on.exit(close(con))
  matrix_data <- read.table(con, header = FALSE, sep = "\t", check.names = FALSE, colClasses = "numeric", na.strings = c("NA", ".", ""))
  R <- as.matrix(matrix_data)
  dimnames(R) <- NULL
  if (!identical(dim(R), rep(nrow(variants), 2L)) || any(!is.finite(R))) stop("LD must be a complete finite square matrix aligned to every variant.")
  if (max(abs(R-t(R))) > qc$matrix_tolerance || max(abs(diag(R)-1)) > qc$matrix_tolerance || max(abs(R)) > 1+qc$matrix_tolerance) stop("LD must be a symmetric correlation matrix with unit diagonal.")
  decomposition <- eigen(R, symmetric = TRUE)
  tolerance <- qc$matrix_tolerance*max(1, max(abs(decomposition$values)))
  if (min(decomposition$values) < -tolerance) stop("LD matrix is not positive semidefinite.")
  # Preserve the supplied R; no shrinkage, allele flips or outlier removal.
  attr(R, "eigen") <- decomposition
  hashes <- c(variants = external_sha256(variant_file), ld = external_sha256(ld_file))
  if (!identical(unname(hashes), tolower(c(entry$variants_sha256, entry$ld_sha256)))) stop("Input files changed while reading.")
  list(variants = variants, R = R, paths = c(variant_file, ld_file),
    hashes = hashes,
    window = window, minimum_eigenvalue = min(decomposition$values), eigenvalue_tolerance = tolerance)
}

external_diagnostics <- function(variants, R, n, qc = EXTERNAL_QC) {
  z <- variants$beta/variants$se
  if (any(!is.finite(z))) stop("Nonfinite beta/SE z score.")
  s <- susieR::estimate_s_rss(z = z, R = R, n = n)
  if (length(s) != 1L || !is.finite(s) || s < 0 || s > 1) stop("Nonfinite or invalid estimate_s_rss diagnostic.")
  kriging <- susieR::kriging_rss(z = z, R = R, n = n, s = s)
  d <- kriging$conditional_dist
  fields <- c("z", "condmean", "condvar", "z_std_diff", "logLR")
  if (!is.data.frame(d) || nrow(d) != length(z) || !all(fields %in% names(d)) || any(!is.finite(as.matrix(d[, fields]))) || any(d$condvar <= 0)) stop("Invalid kriging diagnostics; inference unavailable.")
  d <- data.frame(variant_id = variants$variant_id, variant_type = variants$variant_type, raw_z = z, d,
    severe_outlier = abs(d$z_std_diff) > qc$kriging_abs_residual_warning | (d$logLR > qc$kriging_logLR_warning & abs(d$z) > qc$kriging_abs_z_for_logLR), stringsAsFactors = FALSE)
  list(s = s, table = d, warning = s > qc$s_warning || any(d$severe_outlier),
    outliers = sum(d$severe_outlier), max_abs_standardized_residual = max(abs(d$z_std_diff)))
}
external_fit <- function(z, R, n, qc = EXTERNAL_QC) {
  susieR::susie_rss(z = z, R = R, n = n, estimate_residual_variance = FALSE,
    L = min(10L, length(z)), max_iter = qc$max_iter, coverage = qc$coverage,
    n_purity = length(z))
}
external_usable_fit <- function(fit, p) {
  if (!isTRUE(fit$converged)) stop("SuSiE did not converge; PIPs and credible sets are unavailable.")
  pip <- susieR::susie_get_pip(fit)
  if (length(pip) != p || any(!is.finite(pip) | pip < 0 | pip > 1)) stop("Invalid SuSiE posterior inclusion probabilities.")
  pip
}
external_empty_cs <- function() data.frame(cs_id = character(), component = integer(), n_variants = integer(), posterior_mass = numeric(),
  min_abs_corr = numeric(), mean_abs_corr = numeric(), median_abs_corr = numeric(), first_pos = numeric(), last_pos = numeric(),
  boundary_sensitive = logical(), status = character(), variant_ids = character(), stringsAsFactors = FALSE)
external_credible_sets <- function(fit, variants, R, window, status, qc = EXTERNAL_QC) {
  cs <- susieR::susie_get_cs(fit, Xcorr = R, coverage = qc$coverage, n_purity = nrow(variants))
  output <- external_empty_cs()
  membership <- rep("", nrow(variants))
  for (i in seq_along(cs$cs)) {
    members <- cs$cs[[i]]; name <- names(cs$cs)[i]
    if (is.na(name) || !nzchar(name)) name <- paste0("CS", i)
    edge <- variants$pos[members]-window[1] <= qc$boundary_distance_bp | window[2]-variants$pos[members] <= qc$boundary_distance_bp
    purity <- cs$purity[i, , drop = FALSE]
    output <- rbind(output, data.frame(cs_id = name, component = cs$cs_index[i], n_variants = length(members),
      posterior_mass = cs$coverage[i], min_abs_corr = purity$min.abs.corr, mean_abs_corr = purity$mean.abs.corr,
      median_abs_corr = purity$median.abs.corr, first_pos = min(variants$pos[members]), last_pos = max(variants$pos[members]),
      boundary_sensitive = any(edge), status = status, variant_ids = paste(variants$variant_id[members], collapse = ";")))
    membership[members] <- ifelse(nzchar(membership[members]), paste(membership[members], name, sep = ";"), name)
  }
  list(table = output, membership = membership)
}

external_plot <- function(variants, pip, summary, path) {
  png(path, width = 1800, height = 1150, res = 150, type = "cairo")
  on.exit(dev.off())
  par(mfrow = c(2, 1), mar = c(4, 5.2, 3.5, 1.2), oma = c(3.1, 0, 4, 0), family = "sans")
  palette <- c(SNP = "#22679a", insertion = "#c87512", deletion = "#813d9c", complex_indel = "#237843")
  symbols <- c(SNP = 16, insertion = 17, deletion = 15, complex_indel = 18)
  if (is.null(variants) || !nrow(variants)) {
    for (i in 1:2) { plot.new(); text(.5, .5, if (i == 1) "External estimates unavailable: input validation failed" else "No model PIPs", cex = 1.1) }
  } else {
    x <- variants$pos/1e6; low <- variants$beta-qnorm(.975)*variants$se; high <- variants$beta+qnorm(.975)*variants$se
    colors <- unname(palette[variants$variant_type]); pch <- unname(symbols[variants$variant_type])
    plot(x, variants$beta, ylim = range(c(0, low, high)), type = "n", xlab = paste0(variants$chrom[1], " position (Mb; GRCh38)"), ylab = "External GWAS ALT log-odds ratio", main = "External GWAS estimates and pointwise normal-approximation 95% intervals")
    abline(h = 0, col = "#bbc3cb", lty = 2)
    segments(x, low, x, high, col = adjustcolor(colors, alpha.f = .35)); points(x, variants$beta, pch = pch, col = colors, cex = .6)
    present <- names(palette)[names(palette) %in% variants$variant_type]
    legend("topright", legend = present, col = palette[present], pch = symbols[present], bty = "n", horiz = TRUE)
    plot(x, rep(0, length(x)), ylim = c(0, 1), type = "n", xlab = paste0(variants$chrom[1], " position (Mb; GRCh38)"), ylab = "Exploratory posterior inclusion probability", main = sprintf("SuSiE-RSS with independent reference LD; effective N = %.2f", summary$primary_n))
    if (any(is.finite(pip))) points(x, pip, pch = pch, col = colors, cex = .8)
    else text(mean(range(x)), .5, "PIPs unavailable; see status and diagnostics", col = "#864316")
  }
  mtext(paste(summary$locus_id, "| WINDOW-LIMITED EXTERNAL FINE-MAPPING"), outer = TRUE, line = 2.3, font = 2, cex = 1.15)
  mtext(summary$status, outer = TRUE, line = .8, col = "#864316", cex = .9)
  mtext("Exploratory model values, not validated causal variants or cohort haplotype effects.", side = 1, outer = TRUE, line = .7, cex = .8)
  mtext("Binary-trait RSS approximation; independent LD may differ from the GWAS. Window boundaries can truncate signals.", side = 1, outer = TRUE, line = 1.9, cex = .75)
}

external_run_locus <- function(input_dir, entry, manifest, out, qc = EXTERNAL_QC) {
  folder <- file.path(out, entry$locus_id); dir.create(folder)
  summary <- list(locus_id = entry$locus_id, status = "unavailable_input_validation", reason = NULL,
    analysis_scope = "WINDOW-LIMITED_EXTERNAL", validated_causal_candidates = FALSE,
    total_n = manifest$total_n, effective_n = manifest$effective_n, binary_trait_RSS_approximation = TRUE,
    sample_size_basis = "Nominal study-level case/control totals; no per-variant measured sample sizes were supplied.",
    primary_n = manifest$effective_n, primary_n_basis = "4 / (1/n_cases + 1/n_controls)", sensitivity_n = manifest$total_n,
    window_start = entry$window_start, window_end = entry$window_end, n_variants = 0L,
    n_snps = 0L, n_insertions = 0L, n_deletions = 0L, n_complex_indels = 0L, converged = FALSE,
    credible_sets = 0L, boundary_sensitive_credible_sets = 0L, s = NA_real_, diagnostic_outliers = NA_integer_,
    sensitivity_status = "not_run", sensitivity_max_abs_pip_difference = NA_real_, sensitivity_credible_sets = NA_integer_, warnings = character())
  loaded <- external_capture(external_load_locus(input_dir, entry, manifest$total_n, qc))
  variants <- NULL; pip <- numeric(); sensitivity_pip <- numeric(); cs <- external_empty_cs(); sensitivity_cs <- external_empty_cs(); membership <- character()
  diagnostics <- data.frame(variant_id = character(), variant_type = character(), raw_z = numeric(), z = numeric(), condmean = numeric(), condvar = numeric(), z_std_diff = numeric(), logLR = numeric(), severe_outlier = logical())
  if (!is.null(loaded$error)) summary$reason <- loaded$error else {
    data <- loaded$value; variants <- data$variants; p <- nrow(variants); pip <- rep(NA_real_, p); sensitivity_pip <- pip; membership <- rep(NA_character_, p)
    summary$n_variants <- p; summary$n_snps <- sum(variants$variant_type == "SNP"); summary$n_insertions <- sum(variants$variant_type == "insertion"); summary$n_deletions <- sum(variants$variant_type == "deletion")
    summary$n_complex_indels <- sum(variants$variant_type == "complex_indel")
    summary$minimum_external_p <- min(variants$pval)
    summary$strong_external_association_in_window <- any(variants$pval <= qc$strong_association_p)
    summary$input_sha256 <- as.list(data$hashes); summary$minimum_LD_eigenvalue <- data$minimum_eigenvalue; summary$LD_eigenvalue_tolerance <- data$eigenvalue_tolerance
    summary$strong_boundary_signal <- any(variants$pval <= qc$strong_association_p & (variants$pos-data$window[1] <= qc$boundary_distance_bp | data$window[2]-variants$pos <= qc$boundary_distance_bp))
    diagnostic_result <- external_capture(external_diagnostics(variants, data$R, manifest$effective_n, qc))
    summary$warnings <- c(loaded$warnings, diagnostic_result$warnings)
    if (!is.null(diagnostic_result$error)) { summary$status <- "unavailable_LD_diagnostics"; summary$reason <- diagnostic_result$error } else {
      d <- diagnostic_result$value; diagnostics <- d$table; summary$s <- d$s; summary$diagnostic_outliers <- d$outliers; summary$max_abs_standardized_residual <- d$max_abs_standardized_residual
      summary$status <- "unavailable_fit"
      primary <- external_capture({
        if (d$s >= qc$s_block) stop("Conservative diagnostic block: estimate_s_rss >= ", qc$s_block, "; severe GWAS/LD incompatibility, no PIPs reported. This is a predeclared heuristic, not a universally calibrated threshold.")
        fit <- external_fit(variants$beta/variants$se, data$R, manifest$effective_n, qc)
        values <- external_usable_fit(fit, p)
        status <- if (d$warning) "exploratory_with_LD_warning" else "exploratory_window_limited"
        sets <- external_credible_sets(fit, variants, data$R, data$window, status, qc)
        list(fit = fit, pip = values, sets = sets, status = status)
      })
      summary$warnings <- c(summary$warnings, primary$warnings)
      if (!is.null(primary$error)) {
        summary$reason <- primary$error
        if (d$s >= qc$s_block) summary$status <- "unavailable_severe_LD_incompatibility"
      } else {
        fit <- primary$value$fit; pip <- primary$value$pip; cs <- primary$value$sets$table; membership <- primary$value$sets$membership
        summary$status <- primary$value$status; summary$converged <- TRUE; summary$iterations <- fit$niter
        summary$credible_sets <- nrow(cs); summary$boundary_sensitive_credible_sets <- sum(cs$boundary_sensitive)
        saveRDS(fit, file.path(folder, "model_effective_n.rds"))
        sensitivity <- external_capture({
          alternative_diagnostics <- external_diagnostics(variants, data$R, manifest$total_n, qc)
          if (alternative_diagnostics$s >= qc$s_block) stop("Nominal total-N sensitivity failed the predeclared severe LD incompatibility gate.")
          alternative <- external_fit(variants$beta/variants$se, data$R, manifest$total_n, qc)
          alternative_pip <- external_usable_fit(alternative, p)
          alternative_status <- if (alternative_diagnostics$warning) "exploratory_with_LD_warning" else "exploratory_window_limited"
          alternative_cs <- external_credible_sets(alternative, variants, data$R, data$window, alternative_status, qc)
          list(fit = alternative, pip = alternative_pip, cs = alternative_cs$table, diagnostics = alternative_diagnostics, status = alternative_status)
        })
        summary$warnings <- c(summary$warnings, sensitivity$warnings)
        if (!is.null(sensitivity$error)) { summary$sensitivity_status <- "unavailable_fit"; summary$sensitivity_reason <- sensitivity$error } else {
          sensitivity_pip <- sensitivity$value$pip; sensitivity_cs <- sensitivity$value$cs
          summary$sensitivity_status <- "converged"; summary$sensitivity_max_abs_pip_difference <- max(abs(pip-sensitivity_pip)); summary$sensitivity_credible_sets <- nrow(sensitivity_cs)
          summary$sensitivity_model_status <- sensitivity$value$status; summary$sensitivity_s <- sensitivity$value$diagnostics$s
          summary$sensitivity_diagnostic_outliers <- sensitivity$value$diagnostics$outliers
          external_write_table(sensitivity$value$diagnostics$table, file.path(folder, "diagnostics_total_n.tsv"))
          saveRDS(sensitivity$value$fit, file.path(folder, "model_total_n.rds"))
        }
      }
    }
    if (!identical(unname(vapply(data$paths, external_sha256, character(1))), unname(data$hashes))) {
      summary$status <- "unavailable_input_changed"; summary$reason <- "Input files changed during analysis; all model values discarded."
      summary$converged <- FALSE; summary$credible_sets <- 0L; summary$boundary_sensitive_credible_sets <- 0L
      summary$sensitivity_status <- "discarded_input_changed"; summary$sensitivity_max_abs_pip_difference <- NA_real_; summary$sensitivity_credible_sets <- NA_integer_
      pip[] <- NA_real_; sensitivity_pip[] <- NA_real_; membership[] <- NA_character_; cs <- external_empty_cs(); sensitivity_cs <- external_empty_cs()
      unlink(file.path(folder, c("model_total_n.rds", "model_effective_n.rds")))
    }
  }
  if (is.null(variants)) variants <- as.data.frame(setNames(rep(list(character()), length(EXTERNAL_FIELDS)), EXTERNAL_FIELDS), stringsAsFactors = FALSE)
  output <- cbind(variants, data.frame(pip = pip, pip_total_n = sensitivity_pip, credible_set_ids = membership,
    model_status = rep(summary$status, nrow(variants)), stringsAsFactors = FALSE))
  external_write_table(output, file.path(folder, "variant_pip.tsv")); external_write_table(cs, file.path(folder, "credible_sets.tsv"))
  external_write_table(sensitivity_cs, file.path(folder, "credible_sets_total_n.tsv")); external_write_table(diagnostics, file.path(folder, "diagnostics.tsv"))
  external_plot(if (nrow(variants)) variants else NULL, pip, summary, file.path(folder, "external_finemapping.png"))
  summary$warnings <- unique(summary$warnings)
  external_json(summary, file.path(folder, "summary.json"))
  writeLines(c(paste(entry$locus_id, "WINDOW-LIMITED EXTERNAL FINE-MAPPING"), paste("Status:", summary$status),
    if (!is.null(summary$reason)) paste("Reason:", summary$reason),
    "These external GWAS estimates and model PIPs do not estimate cohort haplotype contrasts or establish causality.",
    if (summary$converged) sprintf("The primary model jointly fitted %d SNPs, %d insertions, %d deletions and %d complex indels within the supplied target window.", summary$n_snps, summary$n_insertions, summary$n_deletions, summary$n_complex_indels) else "No usable primary SuSiE model is reported: input validation, diagnostics or model fitting did not complete successfully; consult the status and reason.",
    "Variants outside the window or absent from the reference panel are not covered.",
    "The independent LD matrix is used without shrinkage, automatic allele flipping, or removal of diagnostic outliers.",
    "PIPs and 95% credible sets are exploratory and conditional on the variants, LD, model and window; no credible set is a validated causal-candidate set.",
    "External beta intervals use beta +/- 1.959964 * SE, are pointwise normal approximations, and are distinct from SuSiE credible sets.",
    "The binary-trait z=beta/SE RSS approximation uses balanced case-control effective N=4/(1/n_cases+1/n_controls), with a separate nominal total-N sensitivity fit. Neither N is a measured per-variant sample count.",
    "Credible sets near the window boundary are flagged; wider-region analysis is required before treating a locus as resolved."), file.path(folder, "INTERPRETATION.txt"))
  summary
}

run_external_finemapping <- function(input_dir, out, workers = 4L, locus = NULL) {
  for (package in c("susieR", "jsonlite", "openssl")) if (!requireNamespace(package, quietly = TRUE)) stop("Required R package missing: ", package)
  if (length(workers) != 1L || !is.finite(workers) || workers < 1 || workers != floor(workers)) stop("workers must be a positive integer.")
  manifest <- external_manifest(input_dir)
  script_hash <- external_sha256(EXTERNAL_SCRIPT_PATH)
  entries <- manifest$data$loci
  if (!is.null(locus)) {
    entries <- Filter(function(x) x$locus_id == locus, entries)
    if (length(entries) != 1L) stop("Requested locus is not uniquely present in the input manifest.")
  }
  if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) stop("Use a new or empty output directory; previous runs are preserved.")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  file.copy(manifest$path, file.path(out, "source_input_manifest.json"), overwrite = FALSE)
  file.copy(EXTERNAL_SCRIPT_PATH, file.path(out, "finemap_external.R"), overwrite = FALSE)
  # Fixed per-locus computations and exact CS purity avoid worker-dependent sampling.
  summaries <- parallel::mclapply(entries, function(entry) external_run_locus(input_dir, entry, manifest, out),
    mc.cores = min(as.integer(workers), length(entries)), mc.preschedule = FALSE)
  if (any(vapply(summaries, inherits, logical(1), "try-error"))) stop("A worker failed outside the scientific status handling; run is incomplete.")
  if (external_sha256(manifest$path) != manifest$sha256) stop("Input manifest changed during analysis; run is incomplete.")
  if (external_sha256(EXTERNAL_SCRIPT_PATH) != script_hash) stop("Fitting script changed during analysis; run is incomplete.")
  columns <- c("locus_id", "status", "n_variants", "n_snps", "n_insertions", "n_deletions", "n_complex_indels", "converged", "credible_sets", "boundary_sensitive_credible_sets", "s", "diagnostic_outliers", "sensitivity_status", "sensitivity_max_abs_pip_difference", "sensitivity_credible_sets")
  table <- do.call(rbind, lapply(summaries, function(x) as.data.frame(x[columns], stringsAsFactors = FALSE)))
  external_write_table(table, file.path(out, "locus_summary.tsv"))
  qc <- EXTERNAL_QC
  record <- list(analysis = "exploratory_window_limited_external_SuSiE_RSS", private_results = TRUE,
    input_dir = normalizePath(input_dir), input_manifest_sha256 = manifest$sha256, source_provenance = manifest$data$provenance,
    script_sha256 = script_hash,
    requested_locus = locus, loci_attempted = length(entries), loci_total = length(manifest$data$loci), workers = as.integer(workers),
    n_cases = manifest$data$n_cases, n_controls = manifest$data$n_controls, total_n = manifest$total_n, effective_n = manifest$effective_n,
    sample_size_basis = "Nominal study-level case/control totals; gwas_n is not a measured per-variant sample size.",
    primary_n = manifest$effective_n, primary_n_basis = "4 / (1/n_cases + 1/n_controls)", sensitivity_n = manifest$total_n,
    sample_size_method_note = "Balanced effective N is the primary binary-trait RSS approximation, also used in published workflows such as eLife 88768; this is not a universal model requirement. Nominal total N is retained as sensitivity analysis.",
    binary_trait_approximation = TRUE, z_score = "ALT beta / SE", residual_variance_estimated = FALSE,
    L = "min(10, number of matched variants)", max_iter = qc$max_iter, coverage = qc$coverage,
    prior_tol = eval(formals(susieR::susie_suff_stat)$prior_tol), LD_shrinkage = FALSE,
    automated_allele_flips = FALSE, diagnostic_outlier_removal = FALSE, thresholds = qc,
    diagnostic_threshold_note = "Conservative reporting heuristics, not calibrated association significance tests or universally validated mismatch thresholds. s>=0.5, invalid/nonfinite diagnostics or unconverged fits have no PIPs; milder finite diagnostic warnings remain explicitly exploratory.",
    validated_causal_candidates = FALSE, cohort_haplotype_inference_modified = FALSE,
    R_version = R.version.string, packages = as.list(vapply(c("susieR", "jsonlite", "openssl"), function(x) as.character(packageVersion(x)), character(1))),
    methods = as.list(EXTERNAL_METHODS), status_counts = as.list(table(table$status)),
    completed_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
  outputs <- list.files(out, recursive = TRUE, full.names = TRUE)
  record$output_sha256 <- as.list(setNames(vapply(outputs, external_sha256, character(1)), substring(outputs, nchar(out)+2L)))
  external_json(record, file.path(out, "run_manifest.json"))
  index <- c("# Window-limited external SNP and indel fine-mapping", "",
    "External GWAS beta estimates and exploratory model probabilities. No cohort haplotype associations or validated causal variants are asserted.", "",
    "| Region | Status | SNPs | Insertions | Deletions | Complex indels | Figure | Model details |", "| --- | --- | ---: | ---: | ---: | ---: | --- | --- |")
  for (x in summaries) index <- c(index, sprintf("| %s | %s | %d | %d | %d | %d | [PNG](%s/external_finemapping.png) | [Interpretation](%s/INTERPRETATION.txt), [PIPs](%s/variant_pip.tsv), [QC](%s/summary.json) |", x$locus_id, x$status, x$n_snps, x$n_insertions, x$n_deletions, x$n_complex_indels, x$locus_id, x$locus_id, x$locus_id, x$locus_id))
  writeLines(index, file.path(out, "index.md")); writeLines("Completed external exploratory fitting; consult per-locus status.", file.path(out, "COMPLETED.txt"))
  invisible(record)
}

if (sys.nframe() == 0L) {
  argv <- commandArgs(TRUE)
  if (length(argv) %% 2L || !length(argv)) stop("Use --input-dir PREPARED --out NEW_DIRECTORY [--workers 4] [--locus LOCUS].")
  keys <- argv[seq(1L, length(argv), 2L)]
  if (anyDuplicated(keys) || any(!keys %in% c("--input-dir", "--out", "--workers", "--locus"))) stop("Unknown or duplicate argument.")
  args <- as.list(argv[seq(2L, length(argv), 2L)]); names(args) <- gsub("-", "_", sub("^--", "", keys))
  if (is.null(args$input_dir) || is.null(args$out)) stop("--input-dir and --out are required.")
  if (!is.null(args$workers)) args$workers <- suppressWarnings(as.numeric(args$workers))
  do.call(run_external_finemapping, args)
}
