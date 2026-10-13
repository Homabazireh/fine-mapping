#!/usr/bin/env Rscript
# Compare saved independent R and Python analyses; outputs aggregate checks only.
entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(entry)), "compare_haplotypes.R"))

canonical_json <- function(x) {
  if (is.list(x) && !length(x)) return(list())
  if (is.list(x)) {
    if (!is.null(names(x))) x <- x[order(names(x))]
    x <- lapply(x, canonical_json)
  }
  x
}

verify_parity <- function(python_dir, r_dir, out, python_audit = NULL, r_audit = NULL) {
  for (folder in c(python_dir, r_dir)) if (!file.exists(file.path(folder, "COMPLETED.txt"))) stop("Incomplete comparison run.")
  summaries <- lapply(c(python_dir, r_dir), function(d) jsonlite::fromJSON(file.path(d, "run_summary.json"), simplifyVector = FALSE))
  for (k in c("group_counts", "locus_status_counts", "selected_loci", "haplotype_carrier_tests_performed", "haplotype_overall_tests_performed")) {
    if (!identical(canonical_json(summaries[[1]][[k]]), canonical_json(summaries[[2]][[k]]))) stop("Summary mismatch: ", k)
  }
  checks <- list()
  check_table <- function(name, left, right) {
    # Python csv.DictWriter quotes JSON cells with standard doubled CSV quotes;
    # the R writer emits literal JSON in tab-delimited fields.
    a <- read.delim(left, header = TRUE, sep = "\t", quote = '"', comment.char = "", check.names = FALSE,
                    colClasses = "character", na.strings = character(), fileEncoding = "UTF-8-BOM")
    b <- read_tsv(right)
    if (!identical(names(a), names(b)) || nrow(a) != nrow(b)) stop("Table dimensions/columns differ: ", name)
    # Overlapping windows are emitted in different, scientifically irrelevant
    # locus order by the two streaming implementations.
    keys <- intersect(c("locus_id", "legacy_locus", "variant_id", "sample", "group", "haplotype_id", "group_a", "group_b", "allele_index"), names(a))
    if (nrow(a) && length(keys)) {
      a <- a[do.call(order, a[keys]), , drop = FALSE]
      b <- b[do.call(order, b[keys]), , drop = FALSE]
      rownames(a) <- rownames(b) <- NULL
    }
    max_error <- 0
    for (k in names(a)) {
      av <- a[[k]]; bv <- b[[k]]
      if (identical(av, bv)) next
      if (k %in% c("lead_p", "external_gwas_p")) {
        ok <- identical(decimal_parts(av), decimal_parts(bv))
      } else if (k %in% c("call_exclusion_reasons", "exclusion_reasons", "allele_index_to_sequence", "groups", "n_resolved", "carriers") &&
                 all(grepl("^[{[]", av)) && all(grepl("^[{[]", bv))) {
        ok <- all(vapply(seq_along(av), function(i) identical(canonical_json(jsonlite::fromJSON(av[i], simplifyVector = FALSE)),
                                                            canonical_json(jsonlite::fromJSON(bv[i], simplifyVector = FALSE))), logical(1)))
      } else if (all(toupper(av) %in% c("TRUE", "FALSE")) && all(toupper(bv) %in% c("TRUE", "FALSE"))) {
        ok <- identical(toupper(av), toupper(bv))
      } else {
        an <- suppressWarnings(as.numeric(av)); bn <- suppressWarnings(as.numeric(bv))
        numeric <- all(!is.na(an) | av == "NA") && all(!is.na(bn) | bv == "NA")
        if (!numeric) ok <- FALSE else {
          tolerance <- 1e-7
          finite <- is.finite(an) & is.finite(bn)
          delta <- abs(an[finite]-bn[finite])
          if (length(delta)) max_error <- max(max_error, delta)
          ok <- identical(is.na(an), is.na(bn)) && identical(is.infinite(an), is.infinite(bn)) &&
            all(sign(an[is.infinite(an)]) == sign(bn[is.infinite(bn)])) &&
            all(delta <= tolerance * pmax(1, abs(an[finite]), abs(bn[finite])))
        }
      }
      if (!isTRUE(ok)) stop("Table mismatch: ", name, "/", k, " (no participant details printed)")
    }
    hash <- tools::md5sum(c(left, right))
    list(table = name, rows = nrow(a), columns = ncol(a), values_match = TRUE,
         raw_bytes_identical = unname(hash[1] == hash[2]), maximum_absolute_numeric_difference = max_error,
         python_md5 = unname(hash[1]), R_md5 = unname(hash[2]))
  }
  files <- c("lead_loci.tsv", "regional_variant_frequencies.tsv.gz", "regional_variant_qc.tsv.gz", "selected_markers.tsv",
             "sample_haplotypes.tsv", "haplotype_group_frequencies.tsv", "haplotype_comparisons.tsv",
             "haplotype_overall_comparisons.tsv", "group_phase_qc.tsv", "selected_marker_genotypes.tsv.gz")
  for (name in files) checks[[length(checks)+1L]] <- check_table(name, file.path(python_dir, name), file.path(r_dir, name))
  # R writes the same aggregate plotting table; Python export may live beside its run.
  py_panel <- file.path(python_dir, "panel_c_data.tsv")
  if (file.exists(py_panel)) checks[[length(checks)+1L]] <- check_table("panel_c_data.tsv", py_panel, file.path(r_dir, "panel_c_data.tsv"))
  if (!is.null(python_audit) && !is.null(r_audit)) for (name in c("external_gwas_leads.tsv", "legacy_locus_audit.tsv", "candidate_variants.tsv"))
    checks[[length(checks)+1L]] <- check_table(name, file.path(python_audit, name), file.path(r_audit, name))
  result <- list(verification = "independent_R_vs_Python", verified_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    synthetic = isTRUE(summaries[[2]]$synthetic), all_checks_passed = TRUE, group_counts = summaries[[2]]$group_counts,
    selected_loci = summaries[[2]]$selected_loci, locus_status_counts = summaries[[2]]$locus_status_counts,
    overall_tests = summaries[[2]]$haplotype_overall_tests_performed, pairwise_tests = summaries[[2]]$haplotype_carrier_tests_performed,
    versions = summaries[[2]]$versions, tables = checks,
    numeric_tolerance = list(relative_and_absolute_general = 1e-7,
      reason = "Python formats many numbers to 8 significant figures; R refines conditional OR/CI roots on log odds to 1e-12."),
    interpretation = "Agreement verifies reproduction of the saved analysis. Untested comparisons remain unavailable; descriptive frequency differences do not establish association or causality.",
    scan_scope_note = "R indexed extraction and Python whole-VCF scanning have different total scan counts; selected-window contents must match.")
  jsonlite::write_json(result, out, pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null", digits = 16)
  cat("PASS: ", length(checks), " saved tables agree; aggregate parity summary: ", out, "\n", sep = "")
  invisible(result)
}

if (sys.nframe() == 0L) {
  argv <- commandArgs(TRUE)
  if ("--help" %in% argv) cat("Rscript code/R/verify_parity.R --python-dir DIR --r-dir DIR --out JSON [--python-audit DIR --r-audit DIR]\n")
  else {
    if (length(argv) %% 2L || !length(argv)) stop("Use --help for arguments.")
    flags <- argv[seq(1, length(argv), 2)]; values <- argv[seq(2, length(argv), 2)]
    if (!all(flags %in% c("--python-dir", "--r-dir", "--out", "--python-audit", "--r-audit"))) stop("Unknown parity option.")
    args <- as.list(values); names(args) <- gsub("-", "_", sub("^--", "", flags))
    do.call(verify_parity, args)
  }
}
