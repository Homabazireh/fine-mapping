#!/usr/bin/env Rscript
# Graphics-only re-render from previously validated, saved regional plot tables.
entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(entry)), "plot_all_regions.R"))

replot_all_regions <- function(run_dir, plot_dir, out, workers = 8L) {
  workers <- suppressWarnings(as.integer(workers))
  if (is.na(workers) || workers < 1L) stop("workers must be a positive integer.")
  if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) stop("Use a new or empty output directory.")
  summary <- jsonlite::fromJSON(file.path(run_dir, "run_summary.json"), simplifyVector = FALSE)
  leads <- read_tsv(file.path(run_dir, "lead_loci.tsv"))
  manifest_path <- file.path(plot_dir, "manifest.json")
  manifest <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
  if (!isTRUE(manifest$inputs_stable_during_plotting) || !isTRUE(manifest$all_clusters_included) ||
      !isTRUE(summary$inputs_stable_during_read)) stop("A completed, validated all-region plot export is required.")
  if (anyDuplicated(leads$locus_id) || nrow(leads) != manifest$selected_loci ||
      nrow(leads) != summary$selected_loci || nrow(leads) != summary$distance_clusters_available)
    stop("Source analysis and saved plot scope disagree.")
  analysis_paths <- file.path(run_dir, names(manifest$input_md5))
  if (!identical(unname(tools::md5sum(analysis_paths)), unname(unlist(manifest$input_md5))))
    stop("Source analysis no longer matches the saved plot export.")
  inputs <- c(manifest_path, file.path(plot_dir, "interpretations.tsv"),
    vapply(manifest$regions, function(r) file.path(plot_dir, r$files$data), character(1)))
  before <- tools::md5sum(inputs)
  if (anyNA(before)) stop("Missing saved regional plotting table.")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  global <- c("all_regional_variant_frequencies.tsv.gz", "interpretations.tsv", "index.md", "index.html")
  if (!all(file.copy(file.path(plot_dir, global), file.path(out, global)))) stop("Could not copy saved aggregate outputs.")
  numeric_fields <- c("pos", "distance_kb", "is_exact_lead_match", "n_group", "n_called", "completeness", "alt_copies", "called_chromosomes", "observed_alt_fraction")
  render <- function(i) {
    region <- manifest$regions[[i]]
    if (!grepl("^[A-Za-z0-9_.-]+$", region$directory)) stop("Invalid saved region directory.")
    lead_index <- match(region$locus_id, leads$locus_id)
    if (is.na(lead_index)) stop("Unknown region in saved export.")
    lead <- as.list(leads[lead_index, ])
    table <- read_tsv(file.path(plot_dir, region$files$data))
    if (!identical(names(table), REGION_PLOT_FIELDS)) stop("Saved plotting schema does not match.")
    for (field in numeric_fields) table[[field]] <- suppressWarnings(as.numeric(table[[field]]))
    rows <- lapply(seq_len(nrow(table)), function(j) as.list(table[j, ]))
    folder <- file.path(out, region$directory); dir.create(folder)
    for (type in c("png", "pdf")) draw_all_region(list(rows = rows), lead, summary, file.path(folder, paste0("region.", type)), type)
    for (field in c("data", "interpretation"))
      if (!file.copy(file.path(plot_dir, region$files[[field]]), file.path(out, region$files[[field]]))) stop("Could not copy saved regional table or interpretation.")
    TRUE
  }
  result <- parallel::mclapply(seq_along(manifest$regions), render, mc.cores = workers, mc.preschedule = FALSE)
  if (!all(vapply(result, isTRUE, logical(1)))) stop("A regional graphics worker failed; output is incomplete.")
  if (!identical(before, tools::md5sum(inputs)) ||
      !identical(unname(tools::md5sum(analysis_paths)), unname(unlist(manifest$input_md5))))
    stop("Saved plotting inputs changed during re-rendering.")
  manifest$created_at_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  manifest$rendered_from_saved_plot_data <- TRUE
  manifest$render_workers <- workers
  manifest$saved_plot_data_unchanged <- TRUE
  manifest$versions <- list(R = R.version.string, jsonlite = as.character(packageVersion("jsonlite")))
  jsonlite::write_json(manifest, file.path(out, "manifest.json"), pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null")
  cat("Re-rendered ", length(manifest$regions), " saved regional figures.\n", sep = "")
  invisible(manifest)
}

if (sys.nframe() == 0L) {
  argv <- commandArgs(TRUE)
  if (length(argv) %% 2L || !length(argv)) stop("Use --run-dir RUN --plot-dir SAVED_PLOTS --out NEW_DIRECTORY [--workers 8].")
  keys <- argv[seq(1L, length(argv), by = 2L)]
  if (anyDuplicated(keys) || !all(keys %in% c("--run-dir", "--plot-dir", "--out", "--workers")) ||
      !all(c("--run-dir", "--plot-dir", "--out") %in% keys)) stop("Invalid re-render arguments.")
  args <- as.list(argv[seq(2L, length(argv), by = 2L)]); names(args) <- gsub("-", "_", sub("^--", "", keys))
  do.call(replot_all_regions, args)
}
