#!/usr/bin/env Rscript
# Native R parallel execution of the same scanner over disjoint region unions.
entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(entry)), "compare_haplotypes.R"))
sequential_scan_vcf <- scan_vcf

partition_windows <- function(loci, workers) {
  # Overlapping windows must stay together: a VCF record is read exactly once.
  ord <- order(vapply(loci, `[[`, character(1), "chrom"), vapply(loci, `[[`, numeric(1), "start"))
  components <- list(); current <- integer(); chrom <- ""; end <- -Inf
  for (i in ord) {
    l <- loci[[i]]
    if (length(current) && (l$chrom != chrom || l$start > end)) {
      components[[length(components)+1L]] <- current; current <- integer()
    }
    current <- c(current, i)
    end <- if (length(current) == 1L) l$end else max(end, l$end)
    chrom <- l$chrom
  }
  if (length(current)) components[[length(components)+1L]] <- current
  workers <- min(workers, length(components))
  buckets <- vector("list", workers); size <- rep(0, workers)
  widths <- vapply(components, function(ix) max(vapply(loci[ix], `[[`, numeric(1), "end")) -
    min(vapply(loci[ix], `[[`, numeric(1), "start")) + 1, numeric(1))
  for (j in order(widths, decreasing = TRUE)) {
    k <- which.min(size); buckets[[k]] <- c(buckets[[k]], components[[j]]); size[k] <- size[k] + widths[j]
  }
  lapply(buckets, sort)
}

merge_scan_tables <- function(paths, out) {
  target <- gzfile(out, "wt"); on.exit(close(target))
  expected_header <- NULL
  for (path in paths) {
    con <- gzfile(path, "rt")
    tryCatch({
      header <- readLines(con, n = 1L, warn = FALSE)
      if (is.null(expected_header)) { expected_header <- header; writeLines(header, target) }
      else if (!identical(header, expected_header)) stop("Worker table header mismatch.")
      repeat { lines <- readLines(con, n = 10000L, warn = FALSE); if (!length(lines)) break; writeLines(lines, target) }
    }, finally = close(con))
  }
}

scan_vcf <- function(args, cohort, loci) {
  workers <- args$workers %||% 1L
  if (workers == 1L) return(sequential_scan_vcf(args, cohort, loci))
  if (.Platform$OS.type == "windows") stop("Parallel scanner requires Unix; use workers=1 on Windows.")
  if (!(file.exists(paste0(args$vcf, ".csi")) || file.exists(paste0(args$vcf, ".tbi"))) || !nzchar(Sys.which("bcftools")))
    stop("Parallel scanner requires an indexed VCF and bcftools; use the sequential implementation otherwise.")
  batches <- partition_windows(loci, workers)
  staging <- file.path(args$out, "scan_workers"); dir.create(staging)
  cat(sprintf("Independent R: %d disjoint-window workers for %d loci\n", length(batches), length(loci)))
  results <- parallel::mclapply(seq_along(batches), function(i) {
    local_args <- args; local_args$out <- file.path(staging, sprintf("worker_%03d", i)); dir.create(local_args$out)
    write_regions(file.path(local_args$out, "regions.bed"), loci[batches[[i]]])
    sink(file.path(local_args$out, "progress.log")); on.exit(sink())
    result <- sequential_scan_vcf(local_args, cohort, loci[batches[[i]]])
    list(loci = result$loci, counts = result$counts, directory = local_args$out)
  }, mc.cores = length(batches), mc.preschedule = FALSE, mc.set.seed = FALSE)
  if (any(vapply(results, inherits, logical(1), "try-error"))) stop("Parallel scanner failed; worker files retained for diagnosis.")
  for (name in c("regional_variant_frequencies.tsv.gz", "regional_variant_qc.tsv.gz"))
    merge_scan_tables(vapply(results, function(r) file.path(r$directory, name), character(1)), file.path(args$out, name))
  merged_loci <- loci
  for (i in seq_along(batches)) merged_loci[batches[[i]]] <- results[[i]]$loci
  fields <- unique(unlist(lapply(results, function(r) names(r$counts))))
  counts <- setNames(lapply(fields, function(k) sum(vapply(results, function(r) r$counts[[k]] %||% 0, numeric(1)))), fields)
  unlink(staging, recursive = TRUE)
  list(loci = merged_loci, counts = counts, scan_scope = "indexed_disjoint_window_unions_parallel_R")
}

if (sys.nframe() == 0L) {
  argv <- commandArgs(TRUE); workers <- 8L
  ix <- which(argv == "--workers")
  if (length(ix)) {
    if (length(ix) != 1L || ix == length(argv)) stop("Supply --workers once with a positive integer.")
    workers <- suppressWarnings(as.numeric(argv[ix+1L])); argv <- argv[-c(ix, ix+1L)]
  }
  if (length(workers) != 1L || !is.finite(workers) || workers < 1 || workers != floor(workers) || workers > 64)
    stop("workers must be an integer from 1 to 64.")
  if ("--help" %in% argv) cat("Parallel R wrapper: --workers N (default 8; maximum 64). Other options follow:\n")
  args <- parse_args(argv)
  if (!is.null(args)) { args$workers <- as.integer(workers); run_analysis(args) }
}
