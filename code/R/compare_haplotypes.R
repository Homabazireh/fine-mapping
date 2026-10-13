#!/usr/bin/env Rscript
# Independent R analysis. bcftools is used only to extract indexed VCF windows;
# genotype filtering, phasing, frequencies, exact tests and graphics run in R.

`%||%` <- function(x, y) if (is.null(x)) y else x

project_root <- function() {
  script <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(script)) {
    path <- normalizePath(sub("^--file=", "", script[1]), mustWork = FALSE)
    root <- dirname(dirname(dirname(path)))
    if (file.exists(file.path(root, "config/settings.json"))) return(root)
  }
  path <- normalizePath(getwd())
  while (!file.exists(file.path(path, "config/settings.json"))) {
    parent <- dirname(path)
    if (parent == path) stop("Run from the project root or use the full script path.")
    path <- parent
  }
  path
}

read_tsv <- function(path) {
  x <- read.delim(path, header = TRUE, sep = "\t", quote = "", comment.char = "",
                  check.names = FALSE, colClasses = "character", na.strings = character(),
                  strip.white = FALSE, fileEncoding = "UTF-8-BOM")
  if (anyDuplicated(names(x))) stop("Duplicate TSV column names: ", path)
  x
}

json <- function(x) if (is.list(x) && !length(x)) "{}" else as.character(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", na = "null", digits = 16))
json_array <- function(x) as.character(jsonlite::toJSON(unname(x), auto_unbox = FALSE, na = "null", digits = 16))
counts_list <- function(x) as.list(table(x))
counter <- function(env, key, n = 1) env[[key]] <- (env[[key]] %||% 0) + n
frac <- function(a, n) if (n > 0) a/n else NA_real_

write_rows <- function(path, rows, fields) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "wt") else file(path, "wt")
  on.exit(close(con))
  writeLines(paste(fields, collapse = "\t"), con)
  for (row in rows) write_row(con, row, fields)
}

write_row <- function(con, row, fields) {
  values <- vapply(fields, function(k) {
    v <- row[[k]]
    if (is.null(v) || !length(v) || is.na(v)) return("NA")
    if (is.numeric(v)) return(format(v, digits = 16, scientific = FALSE, trim = TRUE))
    as.character(v)
  }, character(1))
  writeLines(paste(values, collapse = "\t"), con)
}

defaults <- function(root = project_root(), config = file.path(root, "config/settings.json")) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install R package jsonlite.")
  cfg <- jsonlite::fromJSON(config, simplifyVector = FALSE)
  args <- c(cfg$paths, cfg$analysis)
  for (k in names(cfg$paths)) {
    if (!grepl("^/", args[[k]])) args[[k]] <- file.path(root, args[[k]])
  }
  args$out <- file.path(root, "results/cohort_R")
  args$config <- normalizePath(config)
  args$genome_build <- NULL
  args$gwas_alleles_forward_confirmed <- FALSE
  args$reference_genotypes_confirmed <- FALSE
  args$figure <- TRUE
  args$synthetic <- FALSE
  args$metadata <- NULL # explicit --metadata enables fresh metadata preparation
  args$id_map <- NULL
  args
}

parse_args <- function(argv = commandArgs(TRUE)) {
  root <- project_root()
  if ("--help" %in% argv) {
    cat("Usage: Rscript code/R/compare_haplotypes.R --genome-build GRCh38 --out NEW_DIRECTORY\n",
        "Inputs: --config FILE --vcf FILE --gwas FILE --samples FILE\n",
        "Optional fresh group preparation: --metadata FILE --id-map VERIFIED_TSV\n",
        "Quality options mirror config/settings.json, e.g. --min-dp 10 --min-gq 20\n",
        "--test-associations requests tests (skipped unless reference provenance confirmed).\n",
        "--reference-genotypes-confirmed and --gwas-alleles-forward-confirmed require evidence.\n",
        "--no-figure skips the single panel C export; --synthetic labels software fixtures.\n", sep = "")
    return(NULL)
  }
  config <- file.path(root, "config/settings.json")
  ix <- match("--config", argv)
  if (!is.na(ix)) config <- argv[ix + 1L]
  if (!grepl("^/", config)) config <- file.path(root, config)
  args <- defaults(root, config)
  boolean <- c("allow_filter_dot", "allow_subset", "gwas_alleles_forward_confirmed",
               "test_associations", "reference_genotypes_confirmed", "figure", "synthetic")
  numeric <- c("lead_gap", "top_loci", "flank", "max_markers", "min_dp", "min_gq",
               "min_call_rate", "min_mac", "max_indel_length", "min_test_samples",
               "min_resolved_fraction", "min_carriers", "expected_per_group")
  valid <- c(names(args), "metadata", "id_map", "genome_build", "expected_per_group",
             paste0(c("chrom", "pos", "p", "ea", "oa"), "_col"))
  i <- 1L
  while (i <= length(argv)) {
    flag <- argv[i]
    if (!startsWith(flag, "--")) stop("Expected option, got: ", flag)
    key <- gsub("-", "_", sub("^--(no-)?", "", flag))
    if (!key %in% valid) stop("Unknown option: ", flag)
    if (key %in% boolean) {
      args[[key]] <- !startsWith(flag, "--no-")
    } else {
      i <- i + 1L
      if (i > length(argv)) stop("Missing value for ", flag)
      args[[key]] <- if (key %in% numeric) suppressWarnings(as.numeric(argv[i])) else argv[i]
    }
    i <- i + 1L
  }
  args
}

validate_options <- function(args) {
  if (is.null(args$genome_build) || !args$genome_build %in% c("GRCh37", "GRCh38"))
    stop("--genome-build must declare the common build of GWAS and VCF: GRCh37 or GRCh38.")
  if (args$max_markers < 2 || any(!is.finite(unlist(args[c("flank", "lead_gap", "top_loci", "min_dp", "min_gq", "min_mac")]))))
    stop("Invalid marker or quality settings.")
  if (any(unlist(args[c("flank", "lead_gap", "top_loci", "min_dp", "min_gq", "min_mac")]) < 0)) stop("Negative analysis setting.")
  if (!is.finite(args$min_call_rate) || args$min_call_rate <= 0 || args$min_call_rate > 1 ||
      !is.finite(args$min_resolved_fraction) || args$min_resolved_fraction < 0 || args$min_resolved_fraction > 1)
    stop("Invalid call rate or resolved-fraction threshold.")
  if (args$max_indel_length < 1 || args$max_indel_length > 49 || args$min_test_samples < 1 || args$min_carriers < 1)
    stop("Invalid indel length or test count threshold.")
  integer_keys <- c("lead_gap", "top_loci", "flank", "max_markers", "min_mac", "max_indel_length", "min_test_samples", "min_carriers")
  for (k in integer_keys) if (!is.finite(args[[k]]) || args[[k]] != floor(args[[k]])) stop(k, " must be an integer.")
  if (!is.null(args$expected_per_group) && (!is.finite(args$expected_per_group) || args$expected_per_group < 1 || args$expected_per_group != floor(args$expected_per_group)))
    stop("expected_per_group must be a positive integer.")
  decimal_log10(args$p_threshold)
  if (!is.null(args$id_map) && is.null(args$metadata)) stop("--id-map requires --metadata.")
}

read_vcf_header <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  contigs <- list()
  reference <- character()
  samples <- character()
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) break
    if (startsWith(line, "##contig=<")) {
      id <- sub(".*[<,]ID=([^,>]+).*", "\\1", line)
      len <- if (grepl("[<,]length=[0-9]+", line, ignore.case = TRUE))
        as.numeric(sub(".*[<,]length=([0-9]+).*", "\\1", line, ignore.case = TRUE)) else NA_real_
      contigs[[id]] <- len
    }
    if (grepl("^##(reference|assembly)=", line)) reference <- c(reference, line)
    if (startsWith(line, "#CHROM\t")) {
      samples <- strsplit(line, "\t", fixed = TRUE)[[1]][-(1:9)]
      break
    }
    if (!startsWith(line, "#")) break
  }
  if (!length(samples) || anyNA(samples) || anyDuplicated(samples)) stop("VCF requires unique sample IDs.")
  if (!length(contigs)) stop("VCF requires contig header entries.")
  list(samples = samples, contigs = contigs, reference = reference)
}

prepare_groups <- function(metadata, vcf_samples, id_map = NULL) {
  x <- read_tsv(metadata)
  names(x) <- trimws(names(x))
  schemas <- list(c("name", "Diabetes_Group"), c("name (Tubingen ID)", "Category"))
  schemas <- Filter(function(s) all(s %in% names(x)), schemas)
  if (!length(schemas)) stop("Metadata must contain name + Diabetes_Group or name (Tubingen ID) + Category.")
  ids <- x[[schemas[[1]][1]]]
  groups <- x[[schemas[[1]][2]]]
  for (s in schemas) if (!identical(x[[s[1]]], ids) || !identical(x[[s[2]]], groups)) stop("Conflicting metadata schemas.")
  if (anyDuplicated(ids) || any(!nzchar(trimws(ids))) || any(!nzchar(trimws(groups)))) stop("Duplicate or empty metadata values.")
  mapping <- setNames(groups, ids)
  resolved <- mapping[intersect(vcf_samples, ids)]
  if (!is.null(id_map)) {
    aliases <- read_tsv(id_map)
    if (!all(c("metadata_sample", "vcf_sample") %in% names(aliases))) stop("Alias TSV requires metadata_sample, vcf_sample.")
    if (anyDuplicated(aliases$metadata_sample) || anyDuplicated(aliases$vcf_sample)) stop("Aliases must be one to one.")
    if ("status" %in% names(aliases) && any(!tolower(trimws(aliases$status)) %in% c("confirmed", "verified"))) stop("Only verified aliases are accepted.")
    for (i in seq_len(nrow(aliases))) {
      source <- aliases$metadata_sample[i]; target <- aliases$vcf_sample[i]
      if (!source %in% ids || !target %in% vcf_samples) stop("Alias references an unknown participant.")
      if (source != target && (source %in% vcf_samples || target %in% names(resolved))) stop("Alias collides with an exact match.")
      resolved[target] <- mapping[source]
    }
  }
  data.frame(sample = vcf_samples, group = unname(resolved[vcf_samples]), stringsAsFactors = FALSE)
}

read_cohort <- function(args) {
  header <- read_vcf_header(args$vcf)
  opposite <- if (args$genome_build == "GRCh38") "grch37|hg19" else "grch38|hg38"
  if (any(grepl(opposite, header$reference, ignore.case = TRUE))) stop("VCF reference conflicts with declared genome build.")
  x <- if (!is.null(args$metadata)) prepare_groups(args$metadata, header$samples, args$id_map) else read_tsv(args$samples)
  if (!all(c("sample", "group") %in% names(x))) stop("Sample TSV must contain sample and group.")
  if (anyNA(x$sample) || anyNA(x$group) || any(x$sample != trimws(x$sample)) || any(x$group != trimws(x$group)) ||
      any(!nzchar(x$sample)) || any(toupper(x$group) %in% c("", ".", "NA", "N/A", "UNRESOLVED", "UNKNOWN", "TODO"))) stop("Unresolved or whitespace sample/group values.")
  if (anyDuplicated(x$sample)) stop("Duplicate sample IDs.")
  excluded <- setdiff(header$samples, x$sample)
  if (length(setdiff(x$sample, header$samples)) || (length(excluded) && !args$allow_subset)) stop("Sample IDs do not match VCF exactly.")
  groups <- unique(x$group)
  x <- x[match(intersect(header$samples, x$sample), x$sample), , drop = FALSE]
  if (length(groups) < 2 || !nrow(x)) stop("At least two nonempty groups are required.")
  if (!is.null(args$expected_per_group) && any(table(x$group) != args$expected_per_group)) stop("Unexpected group counts.")
  list(header = header, samples = x$sample, mapping = setNames(x$group, x$sample), groups = groups,
       group_counts = setNames(vapply(groups, function(g) sum(x$group == g), integer(1)), groups), excluded = excluded)
}

chrom_number <- function(x) {
  text <- sub("^chr", "", x, ignore.case = TRUE)
  n <- suppressWarnings(as.integer(text))
  n[!grepl("^[0-9]+$", text) | is.na(n) | n < 1 | n > 22] <- NA_integer_
  n
}

# Exact significant digits and decimal exponent are used for thresholds/ranking.
# Floating-point logarithms are only for presentation and checks in tests.
decimal_parts <- function(x) {
  rows <- lapply(as.character(x), function(s) {
    s <- trimws(s)
    if (!grepl("^[+]?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?$", s, perl = TRUE)) stop("Invalid nonnegative decimal p-value: ", s)
    parts <- strsplit(tolower(sub("^\\+", "", s)), "e", fixed = TRUE)[[1]]
    exponent <- if (length(parts) == 2L) as.numeric(parts[2]) else 0
    mantissa <- parts[1]
    integer_digits <- nchar(strsplit(mantissa, ".", fixed = TRUE)[[1]][1])
    digits <- gsub(".", "", mantissa, fixed = TRUE)
    first <- regexpr("[1-9]", digits)[1]
    if (first < 0) return(list(exponent = -Inf, digits = ""))
    significant <- sub("0+$", "", substring(digits, first))
    exponent <- exponent + integer_digits - first
    if (!is.finite(exponent) || exponent > 0 || (exponent == 0 && significant != "1")) stop("P-value outside [0,1]: ", s)
    list(exponent = exponent, digits = significant)
  })
  data.frame(exponent = vapply(rows, `[[`, numeric(1), "exponent"),
             digits = vapply(rows, `[[`, character(1), "digits"), stringsAsFactors = FALSE)
}

decimal_log10 <- function(x) {
  p <- decimal_parts(x)
  vapply(seq_len(nrow(p)), function(i) {
    if (!nzchar(p$digits[i])) return(-Inf)
    first <- substr(p$digits[i], 1, 15)
    log10(as.numeric(first)) - nchar(first) + 1 + p$exponent[i]
  }, numeric(1))
}

decimal_at_most <- function(x, threshold) {
  p <- decimal_parts(x); cutoff <- decimal_parts(threshold)
  p$exponent < cutoff$exponent | (p$exponent == cutoff$exponent & p$digits <= cutoff$digits)
}

select_leads <- function(args, contigs) {
  x <- read_tsv(args$gwas)
  aliases <- list(chrom = c("Chromsome", "Chromosome", "chrom", "chr", "CHR", "CHROM", "#CHROM"),
                  pos = c("Position", "pos", "POS", "BP"), p = c("Pval", "P", "p", "p_value", "pvalue", "P_VALUE"),
                  ea = c("EffectAllele", "effect_allele", "EA", "A1"), oa = c("NonEffectAllele", "other_allele", "OA", "A2", "NEA"))
  columns <- lapply(names(aliases), function(k) {
    col <- args[[paste0(k, "_col")]] %||% intersect(aliases[[k]], names(x))[1]
    if (is.na(col) || !col %in% names(x)) stop("Missing GWAS column for ", k)
    col
  }); names(columns) <- names(aliases)
  lead <- data.frame(chrom = trimws(x[[columns$chrom]]), pos = suppressWarnings(as.numeric(x[[columns$pos]])),
                     ea = toupper(trimws(x[[columns$ea]])), oa = toupper(trimws(x[[columns$oa]])),
                     p = x[[columns$p]], source_row = seq_len(nrow(x)) + 1L, stringsAsFactors = FALSE)
  p_parts <- decimal_parts(lead$p)
  lead$p_exponent <- p_parts$exponent; lead$p_digits <- p_parts$digits
  if (anyNA(lead$pos) || any(!is.finite(lead$pos)) || any(lead$pos < 1 | lead$pos != floor(lead$pos))) stop("Invalid GWAS position.")
  lead$chr_num <- chrom_number(lead$chrom)
  lead$z_abs <- 0
  if (all(c("Beta", "SE") %in% names(x))) {
    beta <- suppressWarnings(as.numeric(x$Beta)); se <- suppressWarnings(as.numeric(x$SE))
    good <- is.finite(beta) & is.finite(se) & se > 0
    lead$z_abs[good] <- abs(beta[good]/se[good])
  }
  lead$liftover_strand <- if ("liftover_strand" %in% names(x)) x$liftover_strand else "."
  lead$liftover_strand[lead$liftover_strand == ""] <- "."
  id_col <- intersect(c("SNP", "rsid", "ID"), names(x))[1]
  lead$input_id <- if (!is.na(id_col)) x[[id_col]] else "."
  counts <- list(input_rows = nrow(x), above_p_threshold = 0L, rows_with_contig_name_translation = 0L)
  lead <- lead[!is.na(lead$chr_num), , drop = FALSE]
  for (i in seq_len(nrow(lead))) {
    chrom <- lead$chrom[i]
    if (!chrom %in% names(contigs)) {
      candidate <- names(contigs)[which(chrom_number(names(contigs)) == lead$chr_num[i])]
      if (length(candidate) != 1L) stop("Cannot uniquely map GWAS contig: ", chrom)
      lead$chrom[i] <- chrom <- candidate
      counts$rows_with_contig_name_translation <- counts$rows_with_contig_name_translation + 1L
    }
    if (!is.na(contigs[[chrom]]) && lead$pos[i] > contigs[[chrom]]) stop("GWAS position outside contig.")
  }
  eligible <- decimal_at_most(lead$p, args$p_threshold)
  counts$above_p_threshold <- sum(!eligible)
  lead <- lead[eligible & grepl("^[ACGT]+$", lead$ea) & grepl("^[ACGT]+$", lead$oa) & lead$ea != lead$oa, , drop = FALSE]
  if (anyDuplicated(paste(lead$chrom, lead$pos, lead$ea, lead$oa))) stop("Duplicate significant GWAS allele key.")
  if (any(!lead$liftover_strand %in% c(".", "+", "-"))) stop("Invalid liftover strand.")
  counts$eligible_significant_rows <- nrow(lead)
  if (!nrow(lead)) stop("No eligible significant autosomal GWAS leads.")
  lead <- lead[order(lead$chr_num, lead$pos, lead$source_row), , drop = FALSE]
  cluster <- cumsum(c(TRUE, lead$chrom[-1] != lead$chrom[-nrow(lead)] | diff(lead$pos) > args$lead_gap))
  chosen <- lapply(split(seq_len(nrow(lead)), cluster), function(ix) {
    rows <- lead[ix, , drop = FALSE]
    best <- rows[order(rows$p_exponent, rows$p_digits, -rows$z_abs, rows$chr_num, rows$pos, rows$ea, rows$oa, rows$source_row)[1], , drop = FALSE]
    best$cluster_rows <- nrow(rows); best$first_sig <- rows$pos[1]; best$last_sig <- tail(rows$pos, 1)
    best
  })
  chosen <- do.call(rbind, chosen)
  chosen <- chosen[order(chosen$p_exponent, chosen$p_digits, -chosen$z_abs, chosen$chr_num, chosen$pos, chosen$ea, chosen$oa, chosen$source_row), , drop = FALSE]
  n_clusters <- nrow(chosen)
  if (args$top_loci > 0) chosen <- head(chosen, args$top_loci)
  loci <- lapply(seq_len(nrow(chosen)), function(i) {
    l <- as.list(chosen[i, , drop = FALSE]); l$ident <- sprintf("L%04d", i)
    l$start <- max(1, l$pos - args$flank)
    l$end <- min(l$pos + args$flank, contigs[[l$chrom]], na.rm = TRUE)
    l$status <- "pending"
    if (!args$gwas_alleles_forward_confirmed) {
      if (l$liftover_strand == "-") l$status <- "negative_liftover_allele_orientation_unconfirmed"
      else if (nchar(l$ea) == 1 && nchar(l$oa) == 1 && paste(sort(c(l$ea, l$oa)), collapse = "") %in% c("AT", "CG"))
        l$status <- "palindromic_lead_orientation_unconfirmed"
    }
    l$anchors <- list(); l$nearby <- list(); l$counts <- new.env(parent = emptyenv())
    l
  })
  list(loci = loci, counts = counts, columns = columns, cluster_count = n_clusters)
}

parse_call <- function(text, format, n_alleles, min_dp, min_gq) {
  fields <- strsplit(text, ":", fixed = TRUE)[[1]]
  get <- function(k) { i <- match(k, format); if (is.na(i) || i > length(fields)) "." else fields[i] }
  fail <- function(reason) list(gt = NULL, phased = FALSE, ps = NA_character_, reason = reason)
  gt <- get("GT")
  if (!nzchar(gt) || grepl(".", gt, fixed = TRUE)) return(fail("missing_or_partial_GT"))
  separators <- nchar(gt) - nchar(gsub("[/|]", "", gt))
  if (separators != 1L) return(fail("not_diploid"))
  if (!grepl("^[0-9]+[/|][0-9]+$", gt)) return(fail("invalid_GT"))
  parts <- strsplit(gt, "[/|]")[[1]]
  if (length(parts) != 2L) return(fail("not_diploid"))
  if (any(!grepl("^[0-9]+$", parts))) return(fail("invalid_GT"))
  alleles <- suppressWarnings(as.integer(parts))
  if (anyNA(alleles) || any(alleles >= n_alleles)) return(fail("invalid_allele_index"))
  if (!get("FT") %in% c(".", "", "PASS")) return(fail("genotype_filter_failed"))
  for (k in c("DP", "GQ")) {
    threshold <- if (k == "DP") min_dp else min_gq
    if (threshold <= 0) next
    value <- suppressWarnings(as.numeric(get(k)))
    if (is.na(value) && !get(k) %in% c("NaN", "nan")) return(fail(paste0("missing_", k)))
    if (!is.finite(value) || value < threshold) return(fail(paste0("low_or_invalid_", k)))
  }
  ps <- get("PS")
  if (ps %in% c("", ".")) ps <- NA_character_
  list(gt = alleles, phased = grepl("|", gt, fixed = TRUE), ps = ps, reason = "ok")
}

classify_variant <- function(ref, alts, limit) {
  if (!length(alts) || any(!grepl("^[ACGT]+$", c(ref, alts))) || any(ref == alts)) return(NULL)
  kinds <- vapply(alts, function(alt) {
    if (nchar(ref) == 1 && nchar(alt) == 1) "SNP"
    else if (nchar(ref) != nchar(alt) && max(nchar(ref), nchar(alt)) <= limit + 1) "indel"
    else "invalid"
  }, character(1))
  if ("invalid" %in% kinds) NULL else if (length(unique(kinds)) > 1L) "SNP/indel" else kinds[1]
}

resolve_diplotype <- function(calls) {
  fail <- function(reason) list(pair = NULL, status = reason, ps = NA_character_)
  if (!length(calls) || any(vapply(calls, function(c) is.null(c$gt), logical(1)))) return(fail("missing_or_low_quality"))
  het <- Filter(function(c) c$gt[1] != c$gt[2], calls)
  ps <- NA_character_
  if (length(het) > 1L) {
    if (any(!vapply(het, `[[`, logical(1), "phased"))) return(fail("unphased_heterozygote"))
    if (any(vapply(het, function(c) is.na(c$ps), logical(1)))) return(fail("missing_PS"))
    sets <- unique(vapply(het, `[[`, character(1), "ps"))
    if (length(sets) != 1L) return(fail("different_phase_sets"))
    ps <- sets[1]; status <- "same_PS"
  } else {
    if (length(het) && het[[1]]$phased) ps <- het[[1]]$ps
    status <- if (length(het)) "single_het_resolved" else "homozygous_resolved"
  }
  h1 <- vapply(calls, function(c) c$gt[1], integer(1))
  h2 <- vapply(calls, function(c) c$gt[2], integer(1))
  first_diff <- which(h1 != h2)[1]
  if (!is.na(first_diff) && h1[first_diff] > h2[first_diff]) { temp <- h1; h1 <- h2; h2 <- temp }
  list(pair = c(paste(h1, collapse = "|"), paste(h2, collapse = "|")), status = status, ps = ps)
}

REGIONAL_FIELDS <- strsplit("locus_id variant_id chrom pos ref alt type distance_to_lead_bp is_exact_lead_match group n_group n_called call_rate allele_index allele allele_copies called_chromosomes allele_frequency carriers carrier_frequency", " ")[[1]]
QC_FIELDS <- strsplit("locus_id variant_id type n_called n_samples call_rate nonmajor_allele_count eligible_nearby_marker call_exclusion_reasons", " ")[[1]]
MARKER_FIELDS <- strsplit("locus_id marker_number variant_id chrom pos ref alt type is_lead allele_index_to_sequence lead_effect_allele_index pooled_call_rate nonmajor_allele_count", " ")[[1]]
SAMPLE_FIELDS <- strsplit("locus_id sample group status PS n_markers haplotype_1 haplotype_2 haplotype_1_alleles haplotype_2_alleles", " ")[[1]]
FREQ_FIELDS <- strsplit("locus_id haplotype_id allele_indices allele_sequences carries_lead_effect_allele group n_group n_resolved n_unresolved resolved_fraction haplotype_copies resolved_chromosomes haplotype_frequency carriers carrier_frequency", " ")[[1]]
COMP_FIELDS <- strsplit("locus_id haplotype_id group_a group_b n_resolved_a n_resolved_b carriers_a carriers_b carrier_frequency_a carrier_frequency_b test_status conditional_odds_ratio_a_vs_b CI95_low CI95_high p_fisher_two_sided q_BH_all_haplotype_tests carrier_frequency_difference_a_minus_b carrier_difference_CI95_low carrier_difference_CI95_high q_BH_all_pairwise_tests", " ")[[1]]
OVERALL_FIELDS <- strsplit("locus_id haplotype_id groups n_resolved carriers test_status test_method n_tables p_fisher_overall q_BH_all_overall_tests", " ")[[1]]
GROUP_QC_FIELDS <- strsplit("locus_id group n_group n_resolved n_unresolved resolved_fraction exclusion_reasons", " ")[[1]]
GENOTYPE_FIELDS <- strsplit("locus_id sample group variant_id GT PS call_status", " ")[[1]]

write_regions <- function(path, loci) {
  x <- do.call(rbind, lapply(loci, function(l) data.frame(chrom = l$chrom, start = l$start - 1, end = l$end)))
  x <- x[order(x$chrom, x$start, x$end), , drop = FALSE]
  merged <- list()
  for (i in seq_len(nrow(x))) {
    r <- as.list(x[i, ]); j <- length(merged)
    if (j && merged[[j]]$chrom == r$chrom && r$start <= merged[[j]]$end) merged[[j]]$end <- max(merged[[j]]$end, r$end)
    else merged[[j + 1L]] <- r
  }
  writeLines(vapply(merged, function(r) paste(r$chrom, r$start, r$end, sep = "\t"), character(1)), path)
}

scan_vcf <- function(args, cohort, loci) {
  input <- args$vcf; extraction <- NULL
  indexed <- file.exists(paste0(input, ".csi")) || file.exists(paste0(input, ".tbi"))
  if (indexed && nzchar(Sys.which("bcftools"))) {
    extraction <- tempfile("R_selected_windows_", fileext = ".vcf")
    on.exit(unlink(extraction), add = TRUE)
    status <- system2(Sys.which("bcftools"), c("view", "--no-version", "--regions-overlap", "0", "-R",
               shQuote(file.path(args$out, "regions.bed")), "-Ov", "-o", shQuote(extraction), shQuote(input)))
    if (status != 0L) stop("bcftools indexed region extraction failed.")
    input <- extraction
  }
  con <- gzfile(input, "rt"); on.exit(close(con), add = TRUE)
  rh <- gzfile(file.path(args$out, "regional_variant_frequencies.tsv.gz"), "wt"); on.exit(close(rh), add = TRUE)
  qh <- gzfile(file.path(args$out, "regional_variant_qc.tsv.gz"), "wt"); on.exit(close(qh), add = TRUE)
  writeLines(paste(REGIONAL_FIELDS, collapse = "\t"), rh); writeLines(paste(QC_FIELDS, collapse = "\t"), qh)
  ix <- match(cohort$samples, cohort$header$samples)
  gi <- lapply(cohort$groups, function(g) which(cohort$mapping == g)); names(gi) <- cohort$groups
  counts <- new.env(parent = emptyenv()); seen <- new.env(hash = TRUE, parent = emptyenv())
  visited <- character(); last_chrom <- ""; last_pos <- 0
  repeat {
    lines <- readLines(con, n = 500L, warn = FALSE)
    if (!length(lines)) break
    for (line in lines) {
      if (startsWith(line, "#") || !nzchar(line)) next
      f <- strsplit(line, "\t", fixed = TRUE)[[1]]
      if (length(f) != 9L + length(cohort$header$samples)) stop("VCF sample count mismatch.")
      chrom <- f[1]; pos <- suppressWarnings(as.numeric(f[2]))
      if (!is.finite(pos) || pos < 1 || pos != floor(pos)) stop("Invalid VCF position.")
      if (chrom != last_chrom) {
        if (chrom %in% visited) stop("VCF contig reappears; coordinate sorting is required.")
        visited <- c(visited, chrom); last_chrom <- chrom; last_pos <- 0
      }
      if (pos < last_pos) stop("VCF positions are not sorted.")
      last_pos <- pos; counter(counts, "vcf_records_scanned")
      hits <- which(vapply(loci, function(l) l$chrom == chrom && pos >= l$start && pos <= l$end, logical(1)))
      if (!length(hits)) next
      counter(counts, "records_in_windows")
      for (j in hits) { counter(loci[[j]]$counts, "records_in_window"); if (pos == loci[[j]]$pos) counter(loci[[j]]$counts, "records_at_lead_position") }
      ref <- f[4]; alts <- strsplit(f[5], ",", fixed = TRUE)[[1]]; alleles <- c(ref, alts)
      key <- paste(chrom, pos, ref, f[5], sep = ":")
      if (exists(key, envir = seen, inherits = FALSE)) stop("Duplicate VCF variant record: ", key)
      seen[[key]] <- TRUE
      kind <- classify_variant(ref, alts, args$max_indel_length)
      if (is.null(kind)) { counter(counts, "excluded_not_SNP_or_small_indel"); next }
      if (f[7] != "PASS" && !(args$allow_filter_dot && f[7] == ".")) { counter(counts, "excluded_site_FILTER"); next }
      format <- strsplit(f[9], ":", fixed = TRUE)[[1]]
      if (anyDuplicated(format) || !"GT" %in% format) stop("Invalid VCF FORMAT.")
      calls <- lapply(f[9L + ix], parse_call, format = format, n_alleles = length(alleles), min_dp = args$min_dp, min_gq = args$min_gq)
      reasons <- vapply(calls, `[[`, character(1), "reason")
      if (any(reasons %in% c("invalid_GT", "invalid_allele_index"))) stop("Invalid genotype: ", key)
      good <- vapply(calls, function(c) !is.null(c$gt), logical(1)); called <- sum(good)
      copy_counts <- tabulate(unlist(lapply(calls, `[[`, "gt")) + 1L, nbins = length(alleles))
      minor <- sum(copy_counts) - max(copy_counts); rate <- called/length(cohort$samples)
      eligible <- rate >= args$min_call_rate && minor >= args$min_mac
      marker <- list(chrom = chrom, pos = pos, ref = ref, alts = alts, alleles = alleles, key = key, kind = kind,
                     calls = calls, called = called, minor = minor)
      counter(counts, "passing_SNP_indel_records"); counter(counts, "high_quality_genotypes", called)
      counter(counts, "high_quality_homozygous_reference_genotypes", sum(vapply(calls, function(c) identical(c$gt, c(0L, 0L)), logical(1))))
      by_group <- lapply(gi, function(index) {
        gt <- lapply(calls[index[good[index]]], `[[`, "gt")
        list(n = length(gt), copies = tabulate(unlist(gt) + 1L, nbins = length(alleles)),
             carriers = vapply(seq_along(alleles) - 1L, function(a) sum(vapply(gt, function(g) a %in% g, logical(1))), integer(1)))
      })
      for (j in hits) {
        l <- loci[[j]]; counter(l$counts, "passing_SNP_indel_records")
        exact <- pos == l$pos && all(c(l$ea, l$oa) %in% alleles)
        if (exact) l$anchors[[length(l$anchors) + 1L]] <- marker
        if (eligible && pos != l$pos) {
          l$nearby[[length(l$nearby) + 1L]] <- marker
          ord <- order(vapply(l$nearby, function(m) abs(m$pos - l$pos), numeric(1)), vapply(l$nearby, `[[`, numeric(1), "pos"), vapply(l$nearby, `[[`, character(1), "key"))
          l$nearby <- head(l$nearby[ord], args$max_markers - 1L)
        }
        write_row(qh, list(locus_id = l$ident, variant_id = key, type = kind, n_called = called, n_samples = length(cohort$samples),
          call_rate = rate, nonmajor_allele_count = minor, eligible_nearby_marker = as.integer(eligible), call_exclusion_reasons = json(counts_list(reasons[!good]))), QC_FIELDS)
        for (g in cohort$groups) {
          stat <- by_group[[g]]
          for (a in seq_along(alleles)) write_row(rh, list(locus_id = l$ident, variant_id = key, chrom = chrom, pos = pos, ref = ref,
            alt = f[5], type = kind, distance_to_lead_bp = pos - l$pos, is_exact_lead_match = as.integer(exact), group = g,
            n_group = length(gi[[g]]), n_called = stat$n, call_rate = frac(stat$n, length(gi[[g]])), allele_index = a - 1L,
            allele = alleles[a], allele_copies = stat$copies[a], called_chromosomes = 2 * stat$n, allele_frequency = frac(stat$copies[a], 2 * stat$n),
            carriers = stat$carriers[a], carrier_frequency = frac(stat$carriers[a], stat$n)), REGIONAL_FIELDS)
        }
        loci[[j]] <- l
      }
    }
    if ((counts$records_in_windows %||% 0) > 0) cat(sprintf("\rR: scanned %d records in selected windows", counts$records_in_windows))
  }
  cat("\n")
  list(loci = loci, counts = as.list(counts), scan_scope = if (is.null(extraction)) "whole_vcf" else "indexed_selected_windows")
}

wilson <- function(a, n) {
  z <- qnorm(.975); p <- a/n; d <- 1 + z*z/n
  center <- (p + z*z/(2*n))/d
  half <- z * sqrt(p*(1-p)/n + z*z/(4*n*n))/d
  c(center - half, center + half)
}

carrier_test <- function(a, na, b, nb) {
  test <- fisher.test(matrix(c(a, na-a, b, nb-b), nrow = 2, byrow = TRUE), alternative = "two.sided", conf.level = .95)
  # Refine the same conditional noncentral-hypergeometric OR and exact interval
  # on log odds; fisher.test's default odds-scale root tolerance is relatively
  # coarse for sparse tables with very wide intervals.
  total <- a+b; support <- seq.int(max(0, total-nb), min(na, total))
  base <- lchoose(na, support)+lchoose(nb, total-support)
  weights <- function(logodds) { logw <- base+support*logodds; w <- exp(logw-max(logw)); w/sum(w) }
  root <- function(fn) exp(uniroot(fn, c(-700, 700), tol = 1e-12)$root)
  odds <- if (a == min(support)) 0 else if (a == max(support)) Inf else root(function(t) sum(support*weights(t))-a)
  low <- if (a == min(support)) 0 else root(function(t) sum(weights(t)[support >= a])-.025)
  high <- if (a == max(support)) Inf else root(function(t) sum(weights(t)[support <= a])-.025)
  wa <- wilson(a, na); wb <- wilson(b, nb); d <- a/na - b/nb
  list(conditional_odds_ratio_a_vs_b = odds, CI95_low = low, CI95_high = high,
       p_fisher_two_sided = test$p.value, carrier_frequency_difference_a_minus_b = d,
       carrier_difference_CI95_low = d - sqrt((a/na-wa[1])^2 + (wb[2]-b/nb)^2),
       carrier_difference_CI95_high = d + sqrt((wa[2]-a/na)^2 + (b/nb-wb[1])^2))
}

overall_test <- function(carriers, n) {
  # Conditional Freeman-Halton test of a G x 2 carrier table. Fixed row totals
  # and pooled carrier count; enumerate every feasible table, never two copies/person.
  total <- sum(carriers); logden <- lchoose(sum(n), total)
  observed <- sum(lchoose(n, carriers)); p <- 0; nt <- 0L
  enumerate <- function(j, remaining, weight) {
    if (j == length(n)) {
      if (remaining < 0 || remaining > n[j]) return(invisible(NULL))
      value <- weight + lchoose(n[j], remaining)
      nt <<- nt + 1L
      if (value <= observed + 1e-10) p <<- p + exp(value - logden)
      return(invisible(NULL))
    }
    lo <- max(0, remaining - sum(n[(j+1L):length(n)])); hi <- min(n[j], remaining)
    if (lo > hi) return(invisible(NULL))
    for (a in seq.int(lo, hi)) enumerate(j+1L, remaining-a, weight+lchoose(n[j], a))
  }
  enumerate(1L, total, 0)
  list(p = min(1, p), n_tables = nt)
}

test_status <- function(args, n, carriers, totals) {
  if (!args$test_associations) "not_requested"
  else if (!args$reference_genotypes_confirmed) "not_tested_reference_provenance"
  else if (min(n) < args$min_test_samples) "too_few_resolved_samples"
  else if (min(n/totals) < args$min_resolved_fraction) "low_resolved_fraction"
  else if (sum(carriers) == 0 || sum(carriers) == sum(n)) "no_carrier_variation"
  else if (min(sum(carriers), sum(n-carriers)) < args$min_carriers) "too_few_carriers_or_noncarriers"
  else "tested_unadjusted"
}

analyse_loci <- function(args, cohort, loci) {
  rows <- setNames(rep(list(list()), 7), c("markers", "samples", "frequencies", "comparisons", "overall", "group_qc", "genotypes"))
  add <- function(name, value) rows[[name]][[length(rows[[name]]) + 1L]] <<- value
  for (j in seq_along(loci)) {
    l <- loci[[j]]
    if (l$status != "pending") next
    if (!length(l$anchors)) l$status <- "no_unique_PASS_sequence_allele_match_at_lead"
    else if (length(l$anchors) > 1L) l$status <- "multiple_matching_lead_records"
    else if (l$anchors[[1]]$called/length(cohort$samples) < args$min_call_rate) l$status <- "lead_below_min_call_rate"
    else if (length(l$nearby) < 1L) l$status <- "fewer_than_two_eligible_markers"
    if (l$status != "pending") { loci[[j]] <- l; next }
    anchor <- l$anchors[[1]]; panel <- c(list(anchor), l$nearby)
    panel <- panel[order(vapply(panel, `[[`, numeric(1), "pos"), vapply(panel, `[[`, character(1), "key"))]
    anchor_index <- which(vapply(panel, function(m) m$key == anchor$key, logical(1)))
    ea_index <- match(l$ea, anchor$alleles) - 1L
    l$status <- "analysed"; l$counts$selected_markers <- length(panel)
    for (k in seq_along(panel)) {
      m <- panel[[k]]
      add("markers", list(locus_id = l$ident, marker_number = k, variant_id = m$key, chrom = m$chrom, pos = m$pos,
        ref = m$ref, alt = paste(m$alts, collapse = ","), type = m$kind, is_lead = as.integer(k == anchor_index),
        allele_index_to_sequence = json(as.list(setNames(m$alleles, seq_along(m$alleles)-1L))),
        lead_effect_allele_index = if (k == anchor_index) ea_index else NA_integer_,
        pooled_call_rate = m$called/length(cohort$samples), nonmajor_allele_count = m$minor))
    }
    assignment <- lapply(seq_along(cohort$samples), function(i) resolve_diplotype(lapply(panel, function(m) m$calls[[i]])))
    pooled <- table(unlist(lapply(assignment, `[[`, "pair")))
    haps <- names(pooled)
    if (length(haps)) {
      values <- do.call(rbind, lapply(strsplit(haps, "|", fixed = TRUE), as.integer))
      ord <- do.call(order, c(list(-as.numeric(pooled)), as.data.frame(values)))
      haps <- haps[ord]
    }
    ids <- setNames(sprintf(paste0(l$ident, "_H%03d"), seq_along(haps)), haps)
    sequence <- function(hap) {
      indices <- as.integer(strsplit(hap, "|", fixed = TRUE)[[1]])
      paste(vapply(seq_along(panel), function(k) panel[[k]]$alleles[indices[k]+1L], character(1)), collapse = "|")
    }
    for (i in seq_along(cohort$samples)) {
      a <- assignment[[i]]; sample <- cohort$samples[i]; g <- cohort$mapping[[sample]]
      add("samples", list(locus_id = l$ident, sample = sample, group = g, status = a$status, PS = if (is.na(a$ps)) "." else a$ps,
        n_markers = length(panel), haplotype_1 = if (is.null(a$pair)) NA_character_ else ids[[a$pair[1]]],
        haplotype_2 = if (is.null(a$pair)) NA_character_ else ids[[a$pair[2]]],
        haplotype_1_alleles = if (is.null(a$pair)) NA_character_ else sequence(a$pair[1]),
        haplotype_2_alleles = if (is.null(a$pair)) NA_character_ else sequence(a$pair[2])))
      for (m in panel) {
        c <- m$calls[[i]]
        add("genotypes", list(locus_id = l$ident, sample = sample, group = g, variant_id = m$key,
          GT = if (is.null(c$gt)) "." else paste(c$gt, collapse = if (c$phased) "|" else "/"),
          PS = if (is.na(c$ps)) "." else c$ps, call_status = c$reason))
      }
    }
    resolved <- setNames(lapply(cohort$groups, function(g) {
      selected <- assignment[which(cohort$mapping == g)]
      valid <- vapply(selected, function(a) !is.null(a$pair), logical(1)); n <- sum(valid); total <- length(selected)
      add("group_qc", list(locus_id = l$ident, group = g, n_group = total, n_resolved = n, n_unresolved = total-n,
        resolved_fraction = n/total, exclusion_reasons = json(counts_list(vapply(selected[!valid], `[[`, character(1), "status")))))
      lapply(selected[valid], `[[`, "pair")
    }), cohort$groups)
    ns <- vapply(resolved, length, integer(1))
    for (hap in haps) {
      carrier <- setNames(integer(length(cohort$groups)), cohort$groups)
      hi <- as.integer(strsplit(hap, "|", fixed = TRUE)[[1]])
      for (g in cohort$groups) {
        pairs <- resolved[[g]]; copies <- sum(unlist(pairs) == hap); n <- ns[[g]]
        carrier[[g]] <- sum(vapply(pairs, function(pair) hap %in% pair, logical(1)))
        add("frequencies", list(locus_id = l$ident, haplotype_id = ids[[hap]], allele_indices = hap, allele_sequences = sequence(hap),
          carries_lead_effect_allele = as.integer(hi[anchor_index] == ea_index), group = g, n_group = cohort$group_counts[[g]],
          n_resolved = n, n_unresolved = cohort$group_counts[[g]]-n, resolved_fraction = n/cohort$group_counts[[g]],
          haplotype_copies = copies, resolved_chromosomes = 2*n, haplotype_frequency = frac(copies, 2*n),
          carriers = carrier[[g]], carrier_frequency = frac(carrier[[g]], n)))
      }
      status <- test_status(args, ns, carrier, cohort$group_counts)
      if (status == "tested_unadjusted" && !length(cohort$groups) %in% c(2L, 3L)) status <- "unsupported_overall_group_count"
      overall <- list(locus_id = l$ident, haplotype_id = ids[[hap]], groups = json_array(cohort$groups), n_resolved = json_array(ns),
        carriers = json_array(carrier), test_status = status, test_method = "fisher_freeman_halton_exact", n_tables = NA_integer_, p_fisher_overall = NA_real_)
      if (status == "tested_unadjusted") { test <- overall_test(carrier, ns); overall$n_tables <- test$n_tables; overall$p_fisher_overall <- test$p }
      add("overall", overall)
      for (pair in combn(cohort$groups, 2, simplify = FALSE)) {
        ga <- pair[1]; gb <- pair[2]; na <- ns[[ga]]; nb <- ns[[gb]]; a <- carrier[[ga]]; b <- carrier[[gb]]
        status <- test_status(args, c(na, nb), c(a, b), cohort$group_counts[pair])
        row <- list(locus_id = l$ident, haplotype_id = ids[[hap]], group_a = ga, group_b = gb, n_resolved_a = na, n_resolved_b = nb,
          carriers_a = a, carriers_b = b, carrier_frequency_a = frac(a, na), carrier_frequency_b = frac(b, nb),
          carrier_frequency_difference_a_minus_b = NA_real_, test_status = status)
        if (status == "tested_unadjusted") row <- modifyList(row, carrier_test(a, na, b, nb))
        add("comparisons", row)
      }
    }
    loci[[j]] <- l
  }
  for (family in c("comparisons", "overall")) {
    pfield <- if (family == "comparisons") "p_fisher_two_sided" else "p_fisher_overall"
    qfield <- if (family == "comparisons") "q_BH_all_pairwise_tests" else "q_BH_all_overall_tests"
    ps <- vapply(rows[[family]], function(r) r[[pfield]] %||% NA_real_, numeric(1))
    qs <- rep(NA_real_, length(ps)); good <- is.finite(ps)
    qs[good] <- p.adjust(ps[good], method = "BH")
    for (i in seq_along(qs)) {
      rows[[family]][[i]][[qfield]] <- qs[i]
      if (family == "comparisons") rows[[family]][[i]]$q_BH_all_haplotype_tests <- qs[i]
    }
  }
  files <- c(markers = "selected_markers.tsv", samples = "sample_haplotypes.tsv", frequencies = "haplotype_group_frequencies.tsv",
    comparisons = "haplotype_comparisons.tsv", overall = "haplotype_overall_comparisons.tsv", group_qc = "group_phase_qc.tsv", genotypes = "selected_marker_genotypes.tsv.gz")
  fields <- list(markers = MARKER_FIELDS, samples = SAMPLE_FIELDS, frequencies = FREQ_FIELDS, comparisons = COMP_FIELDS,
                 overall = OVERALL_FIELDS, group_qc = GROUP_QC_FIELDS, genotypes = GENOTYPE_FIELDS)
  for (k in names(files)) write_rows(file.path(args$out, files[[k]]), rows[[k]], fields[[k]])
  list(loci = loci, rows = rows, pairwise_tests = sum(vapply(rows$comparisons, function(r) r$test_status == "tested_unadjusted", logical(1))),
       overall_tests = sum(vapply(rows$overall, function(r) r$test_status == "tested_unadjusted", logical(1))))
}

display_group <- function(g) {
  display <- c(T2D = "Diabetes", diabetes = "Diabetes", Diabetes = "Diabetes", Prediabetes = "Prediabetes",
               prediabetes = "Prediabetes", NORM = "Control", Control = "Control", control = "Control")
  out <- unname(display[g]); out[is.na(out)] <- g[is.na(out)]; out
}

panel_c_data <- function(cohort, loci) {
  rows <- list()
  display_order <- match(display_group(cohort$groups), c("Diabetes", "Prediabetes", "Control"))
  groups <- cohort$groups[order(display_order, cohort$groups)]
  for (l in loci) {
    m <- if (length(l$anchors) == 1L) l$anchors[[1]] else NULL
    for (g in c(groups, "ALL")) {
      ix <- if (g == "ALL") seq_along(cohort$samples) else which(cohort$mapping == g)
      called <- if (is.null(m)) NA_integer_ else sum(vapply(m$calls[ix], function(c) !is.null(c$gt), logical(1)))
      rows[[length(rows)+1L]] <- list(locus_id = l$ident, chrom = l$chrom, lead_pos = l$pos, status = l$status,
        group = g, n_called = called, n_group = length(ix), completeness = called/length(ix))
    }
  }
  rows
}

draw_panel_c <- function(args, cohort, loci, data) {
  groups <- unique(vapply(data, `[[`, character(1), "group")); groups <- groups[groups != "ALL"]
  colors <- c(Diabetes = "#d55e00", Prediabetes = "#009e73", Control = "#0072b2")
  color <- colors[display_group(groups)]; color[is.na(color)] <- "#64748b"
  n <- length(loci); pooled <- Filter(function(r) r$group == "ALL", data)
  rates <- vapply(pooled, `[[`, numeric(1), "completeness")
  best <- if (any(is.finite(rates))) which.max(replace(rates, !is.finite(rates), -Inf)) else NA_integer_
  height <- max(8, 3.6 + n*.34)
  render <- function() {
    par(fig = c(.255, .97, .20, .77), mar = c(2.8, 0, 0, 0), family = "sans", bg = "white", fg = "#233445", col.axis = "#233445", col.lab = "#233445")
    plot(NA_real_, xlim = c(-2, 118), ylim = c(n+.5, .5), xlab = "",
         ylab = "", axes = FALSE, xaxs = "i", yaxs = "i", cex.lab = 1)
    for (i in seq_len(n)) rect(-2, i-.48, 118, i+.48, col = if (!is.na(best) && i == best) "#fff0ca" else if (i %% 2L) "#f5f7fa" else "white", border = NA)
    abline(v = seq(0, 100, 20), col = "#dce3ea", lwd = .7)
    abline(v = 100*args$min_call_rate, col = "#a54536", lty = 2, lwd = 1.2)
    symbols <- c(16, 17, 15, 18, 8, 3)
    for (k in seq_along(groups)) {
      rows <- Filter(function(r) r$group == groups[k], data)
      x <- vapply(rows, `[[`, numeric(1), "completeness")*100
      points(x, seq_len(n) + (k-(length(groups)+1)/2)*.19, pch = symbols[(k-1L) %% length(symbols)+1L], col = color[k], cex = .85)
    }
    points(rates*100, seq_len(n), pch = 18, col = "#233445", cex = .85)
    axis(1, at = seq(0, 100, 20), labels = paste0(seq(0, 100, 20), "%"), lwd = 0, lwd.ticks = 1, cex.axis = .85)
    labels <- vapply(loci, function(l) paste(l$ident, paste0(l$chrom, ":", format(l$pos, big.mark = ",", scientific = FALSE, trim = TRUE)), sep = "   "), character(1))
    axis(2, at = seq_len(n), labels = labels, las = 1, lwd = 0, lwd.ticks = 0, cex.axis = .8)
    text(109, seq_len(n), vapply(pooled, function(r) if (is.na(r$n_called)) "NA" else paste0(r$n_called, "/", r$n_group), character(1)), cex = .8)
    text(109, .1, "Pooled calls", xpd = NA, cex = .75)
    note <- if (is.na(best)) "No unique exact lead genotypes were available for pooled completeness estimates." else
      sprintf("Highlighted: %s, %d/%d (%.1f%%). %d/%d leads meet the pooled requirement.", loci[[best]]$ident, pooled[[best]]$n_called,
              length(cohort$samples), 100*rates[best], sum(rates >= args$min_call_rate, na.rm = TRUE), n)
    # Headings, legend and footnotes use page coordinates to avoid collisions
    # when the number of loci or graphics device changes.
    par(fig = c(0, 1, 0, 1), mar = rep(0, 4), new = TRUE)
    plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
    text(.045, .954, "c   Usable lead genotypes across selected regions", adj = 0, font = 2, cex = 1.45, col = "#172f46")
    text(.045, .912, sprintf("%d participants | %d GWAS-selected regions | %s", length(cohort$samples), n, args$genome_build), adj = 0, cex = 1, col = "#526579")
    text(.045, .877, sprintf("Usable genotype: DP >= %g, GQ >= %g. Dashed line: %.0f%% pooled completeness requirement.", args$min_dp, args$min_gq, 100*args$min_call_rate), adj = 0, cex = .85, col = "#526579")
    legend(.58, .837, xjust = .5, yjust = .5, legend = c(paste0(display_group(groups), " (n = ", cohort$group_counts[groups], ")"), "All participants"),
           col = c(color, "#233445"), pch = c(symbols[seq_along(groups)], 18), horiz = TRUE, bty = "n", cex = .82)
    text(.59, .16, "Participants with usable genotypes at the lead variant", cex = .9, col = "#233445")
    text(.045, .09, note, adj = 0, cex = .8, col = "#40566d")
    text(.045, .059, "Dots measure genotype availability; haplotype comparisons require eligible phased panels.", adj = 0, cex = .8, col = "#526579")
    text(.045, .025, paste(if (args$synthetic) "SYNTHETIC FIXTURE" else "Real cohort", "| Independent R analysis | Layout inspired by Figueroa et al., Nature Genetics (2024), Fig. 1"), adj = 0, cex = .65, col = "#63758a")
  }
  for (type in c("png", "pdf", "svg")) {
    path <- file.path(args$out, paste0("panel_c.", type))
    if (type == "png") png(path, width = 13, height = height, units = "in", res = 600, type = "cairo")
    else if (type == "pdf") pdf(path, width = 13, height = height, useDingbats = FALSE)
    else svg(path, width = 13, height = height)
    tryCatch(render(), finally = dev.off())
  }
}

write_lead_loci <- function(args, loci) {
  rows <- lapply(loci, function(l) list(locus_id = l$ident, chrom = l$chrom, lead_pos = l$pos, lead_effect_allele = l$ea,
    lead_other_allele = l$oa, lead_p = l$p, gwas_source_line = l$source_row, input_variant_id = l$input_id,
    start = l$start, end = l$end, cluster_significant_rows = l$cluster_rows, cluster_first_significant_pos = l$first_sig,
    cluster_last_significant_pos = l$last_sig, status = l$status, records_in_window = l$counts$records_in_window %||% 0,
    records_at_lead_position = l$counts$records_at_lead_position %||% 0, exact_passing_lead_matches = length(l$anchors),
    selected_markers = l$counts$selected_markers %||% 0))
  write_rows(file.path(args$out, "lead_loci.tsv"), rows, names(rows[[1]]))
}

write_evidence_assessment <- function(args, loci, results) {
  # Native R evidence eligibility audit, based only on this independent run.
  # An external GWAS lead alone is not a fine-mapped or causal variant.
  rows <- lapply(loci, function(l) {
    pair <- Filter(function(r) r$locus_id == l$ident && identical(r$test_status, "tested_unadjusted"), results$rows$comparisons)
    sig <- sum(vapply(pair, function(r) !is.na(r$q_BH_all_pairwise_tests) && r$q_BH_all_pairwise_tests < .05, logical(1)))
    list(locus_id = l$ident, chrom = l$chrom, lead_pos = l$pos, effect_allele = l$ea, other_allele = l$oa,
      external_gwas_p = l$p, panel_status = l$status, cohort_pairwise_tests = length(pair), cohort_q_below_0_05 = sig,
      fine_mapping_status = "not_performed_missing_validated_dense_summary_statistics_and_matched_LD",
      causal_status = "not_established", candidate_status = "external_GWAS_lead_only",
      reason = if (l$status == "analysed") "Haplotype association cannot establish a causal variant; independent fine-mapping and functional evidence required."
      else paste("Cohort haplotype comparison unavailable:", l$status))
  })
  write_rows(file.path(args$out, "variant_evidence_assessment.tsv"), rows, names(rows[[1]]))
}

run_analysis <- function(args) {
  validate_options(args)
  cohort <- read_cohort(args)
  selection <- select_leads(args, cohort$header$contigs)
  if (dir.exists(args$out) && length(list.files(args$out, all.files = TRUE, no.. = TRUE))) stop("Output directory is not empty; use a new directory.")
  input_paths <- unique(c(args$vcf, args$gwas, if (is.null(args$metadata)) args$samples else c(args$metadata, args$id_map)))
  before <- file.info(input_paths)[, c("size", "mtime"), drop = FALSE]
  dir.create(args$out, recursive = TRUE, showWarnings = FALSE)
  samples <- lapply(cohort$samples, function(s) list(sample = s, group = cohort$mapping[[s]]))
  write_rows(file.path(args$out, "samples_used.tsv"), samples, c("sample", "group"))
  write_regions(file.path(args$out, "regions.bed"), selection$loci)
  scan <- scan_vcf(args, cohort, selection$loci)
  results <- analyse_loci(args, cohort, scan$loci)
  loci <- results$loci
  write_lead_loci(args, loci)
  panel <- panel_c_data(cohort, loci)
  write_rows(file.path(args$out, "panel_c_data.tsv"), panel, c("locus_id", "chrom", "lead_pos", "status", "group", "n_called", "n_group", "completeness"))
  if (args$figure) draw_panel_c(args, cohort, loci, panel)
  write_evidence_assessment(args, loci, results)
  stable <- identical(before, file.info(input_paths)[, c("size", "mtime"), drop = FALSE])
  warnings <- c("Frequencies condition on observed high-quality resolved calls and can be biased by missingness.",
    "Physical-distance GWAS clusters and exploratory haplotype carrier tests do not establish causal variants.",
    "Participant-level carrier tests are unadjusted for age, sex, ancestry, relatedness and technical batch.",
    "Overall and pairwise BH families are separate across all loci; chromosomes are not independent observations.",
    "Missing genotypes remain missing; marker panels are phenotype blind and phase sets are never stitched.")
  if (!args$reference_genotypes_confirmed) warnings <- c(warnings, "Reference-call provenance unconfirmed: all inferential tests remain unavailable.")
  if ((scan$counts$high_quality_homozygous_reference_genotypes %||% 0) == 0) warnings <- c(warnings, "No high-quality homozygous reference calls in selected windows; investigate variant-only merging/reference-call absence.")
  if (!any(vapply(loci, function(l) l$status == "analysed", logical(1)))) warnings <- c(warnings, "No eligible haplotype panels: group differences and causal variants cannot be inferred from this run. This does not show absence of effects.")
  summary <- list(implementation = "independent_R", version = "1.0.0", completed_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    synthetic = args$synthetic, genome_build_user_declared = args$genome_build, parameters = args,
    reference_header = cohort$header$reference, group_counts = as.list(cohort$group_counts),
    excluded_unlisted_vcf_samples = unname(cohort$excluded), gwas_column_mapping = selection$columns,
    gwas_counts = selection$counts, distance_clusters_available = selection$cluster_count, selected_loci = length(loci),
    locus_status_counts = counts_list(vapply(loci, `[[`, character(1), "status")), scan_scope = scan$scan_scope, scan_counts = scan$counts,
    haplotype_carrier_tests_performed = results$pairwise_tests, haplotype_overall_tests_performed = results$overall_tests,
    inputs_stable_during_read = stable, warnings = warnings, input_files = lapply(seq_along(input_paths), function(i)
      list(path = normalizePath(input_paths[i]), bytes = before$size[i], mtime = format(before$mtime[i], tz = "UTC", usetz = TRUE))),
    small_input_md5 = as.list(tools::md5sum(setdiff(input_paths, args$vcf))),
    versions = list(R = R.version.string, jsonlite = as.character(packageVersion("jsonlite")),
                    bcftools = if (nzchar(Sys.which("bcftools"))) system2(Sys.which("bcftools"), "--version", stdout = TRUE)[1] else NULL))
  jsonlite::write_json(summary, file.path(args$out, "run_summary.json"), pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null", digits = 16)
  writeLines(warnings, file.path(args$out, "INTERPRETATION.txt"))
  capture.output(sessionInfo(), file = file.path(args$out, "sessionInfo.txt"))
  if (!stable) stop("An input changed during the run; discard outputs and rerun.")
  writeLines("Completed independent R analysis. Read run_summary.json and skipped locus statuses.", file.path(args$out, "COMPLETED.txt"))
  cat("Groups: ", json(as.list(cohort$group_counts)), "\nStatuses: ", json(summary$locus_status_counts),
      "\nOverall tests: ", results$overall_tests, "; pairwise tests: ", results$pairwise_tests, "\nOutputs: ", normalizePath(args$out), "\n", sep = "")
  invisible(summary)
}

if (sys.nframe() == 0L) {
  args <- parse_args()
  if (!is.null(args)) run_analysis(args)
}
