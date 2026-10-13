#!/usr/bin/env Rscript
# RStudio: source("code/R/run_configured_analysis.R")
# run_configured_analysis("config/R_workflow.example.json", dry_run = TRUE)
# Change dry_run to FALSE after editing the selected action's paths/settings.
# source() defines helpers only: it never starts an analysis.
R_WORKFLOW_DIRECTORY <- local({
  sourced <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
  path <- if (length(sourced)) tail(sourced, 1L)[[1L]] else {
    sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  }
  dirname(normalizePath(path, mustWork = TRUE))
})

workflow_object <- function(x, allowed, label) {
  if (!is.list(x) || (length(x) && (is.null(names(x)) || any(!nzchar(names(x))) || anyDuplicated(names(x)))))
    stop(label, " must be a JSON object with unique setting names.")
  unknown <- setdiff(names(x), allowed)
  if (length(unknown)) stop("Unknown ", label, " setting: ", paste(unknown, collapse = ", "))
  x
}
workflow_string <- function(x, label) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x) || grepl("[[:cntrl:]]", x))
    stop(label, " must be a nonempty string without control characters.")
  x
}
workflow_path <- function(x, base, label) {
  x <- path.expand(workflow_string(x, label))
  if (!grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", x)) x <- file.path(base, x)
  normalizePath(x, mustWork = FALSE)
}
workflow_options <- function(options, logicals = character(), numbers = character(), strings = character()) {
  options <- workflow_object(options, c(logicals, numbers, strings), "options")
  for (name in names(options)) {
    value <- options[[name]]
    if (name %in% logicals && (!is.logical(value) || length(value) != 1L || is.na(value))) stop(name, " must be true or false.")
    if (name %in% numbers && (!is.numeric(value) || length(value) != 1L || !is.finite(value))) stop(name, " must be a finite number.")
    if (name %in% strings) workflow_string(value, name)
  }
  options
}
workflow_flags <- function(options, logicals = character()) {
  unlist(lapply(names(options), function(name) {
    flag <- gsub("_", "-", name, fixed = TRUE)
    if (name %in% logicals) return(paste0(if (options[[name]]) "--" else "--no-", flag))
    c(paste0("--", flag), if (is.numeric(options[[name]])) format(options[[name]], scientific = FALSE, trim = TRUE, digits = 16) else options[[name]])
  }), use.names = FALSE)
}

configured_analysis_command <- function(config, action = NULL) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install R package jsonlite.")
  config <- normalizePath(config, mustWork = TRUE)
  settings <- workflow_object(jsonlite::fromJSON(config, simplifyVector = FALSE), c("action", "actions"), "workflow")
  available <- c("external_finemapping", "cohort_analysis", "regional_plots", "haplotype_contrasts", "descriptive_haplotypes")
  if (is.null(action)) action <- settings$action
  action <- workflow_string(action, "action")
  if (!action %in% available) stop("Unknown action: ", action)
  actions <- workflow_object(settings$actions, available, "actions")
  entry <- workflow_object(actions[[action]], c("input", "out", "options"), action)
  options <- if (is.null(entry$options)) list() else entry$options
  base <- dirname(config)
  args <- character()
  add_path <- function(name, value) c(paste0("--", gsub("_", "-", name)), workflow_path(value, base, name))
  if (action == "external_finemapping") {
    input <- workflow_object(entry$input, "config", "external input")
    options <- workflow_options(options,
      numbers = c("workers", "max_effects", "coverage", "max_iter", "plot_width", "plot_height", "plot_dpi"),
      strings = c("locus", "variant_types", "plot_formats"))
    args <- add_path("config", input$config)
    # Settings inside this file use its own directory as their path base.
    if (!file.exists(args[2])) stop("External configuration does not exist: ", args[2])
    if (!is.null(entry$out)) args <- c(args, add_path("out", entry$out))
    args <- c(args, workflow_flags(options))
    script <- "finemap_external.R"
  } else if (action == "cohort_analysis") {
    input <- workflow_object(entry$input, c("config", "vcf", "gwas", "samples", "metadata", "id_map"), "cohort input")
    if (is.null(input$config)) stop("cohort_analysis requires input.config (the authoritative settings file).")
    logicals <- c("allow_filter_dot", "allow_subset", "test_associations", "figure", "synthetic")
    options <- workflow_options(options, logicals,
      numbers = c("workers", "lead_gap", "top_loci", "flank", "max_markers", "min_dp", "min_gq", "min_call_rate",
        "min_mac", "max_indel_length", "min_test_samples", "min_resolved_fraction", "min_carriers", "expected_per_group"),
      strings = c("genome_build", "p_threshold", "chrom_col", "pos_col", "p_col", "ea_col", "oa_col"))
    if (is.null(options$genome_build)) stop("cohort_analysis requires options.genome_build.")
    workers <- if (is.null(options$workers)) 1 else options$workers
    if (workers < 1 || workers > 64 || workers != floor(workers)) stop("workers must be an integer from 1 to 64.")
    script <- if (workers == 1) "compare_haplotypes.R" else "compare_haplotypes_parallel.R"
    if (workers == 1) options$workers <- NULL
    args <- unlist(lapply(names(input), function(name) add_path(name, input[[name]])), use.names = FALSE)
    args <- c(args, add_path("out", entry$out), workflow_flags(options, logicals),
      "--no-reference-genotypes-confirmed", "--no-gwas-alleles-forward-confirmed")
  } else {
    input <- workflow_object(entry$input, "run_dir", "plot input")
    args <- c(add_path("run_dir", input$run_dir), add_path("out", entry$out))
    if (action == "regional_plots") {
      options <- workflow_options(options, logicals = "include_indels")
      if (isTRUE(options$include_indels)) args <- c(args, "--include-indels")
      script <- "plot_all_regions.R"
    } else if (action == "descriptive_haplotypes") {
      options <- workflow_options(options, numbers = "columns")
      if (!is.null(options$columns) && (options$columns < 1 || options$columns != floor(options$columns))) stop("columns must be a positive integer.")
      args <- c(args, workflow_flags(options))
      script <- "plot_descriptive_haplotypes.R"
    } else {
      workflow_options(options)
      script <- "plot_haplotype_contrasts.R"
    }
  }
  executable <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  list(action = action, config = config, executable = executable,
    script = file.path(R_WORKFLOW_DIRECTORY, script), args = args)
}

workflow_invoke <- function(command) {
  # system2 passes args through a shell on Unix; quote EACH argument, including
  # script/input/output paths, rather than concatenating unquoted user text.
  status <- system2(command$executable, shQuote(c(command$script, command$args)))
  if (!identical(as.integer(status), 0L)) stop("R analysis failed with exit status ", status, ". Read the child-script error; outputs may be incomplete.")
  invisible(command)
}
run_configured_analysis <- function(config, action = NULL, dry_run = FALSE) {
  if (!is.logical(dry_run) || length(dry_run) != 1L || is.na(dry_run)) stop("dry_run must be TRUE or FALSE.")
  command <- configured_analysis_command(config, action)
  cat("Action: ", command$action, "\n",
    paste(shQuote(c(command$executable, command$script, command$args)), collapse = " "), "\n", sep = "")
  if (dry_run) return(invisible(command))
  workflow_invoke(command)
}
configured_analysis_help <- function() cat(
  "Usage: Rscript code/R/run_configured_analysis.R --config FILE [--action ACTION] [--dry-run]\n",
  "RStudio: source('code/R/run_configured_analysis.R'); run_configured_analysis('config/R_workflow.example.json', dry_run=TRUE)\n",
  "Actions: external_finemapping, cohort_analysis, regional_plots, haplotype_contrasts, descriptive_haplotypes.\n",
  "Edit actions.<selected action>.input/out/options; only the selected action is validated or executed.\n",
  "Workflow paths are relative to FILE. External input.config delegates to its own config (max_effects, coverage, workers, plots).\n",
  "Cohort input.config retains existing settings.json thresholds; its native paths remain relative to the repository.\n",
  "Dry run prints arguments without requiring restored private inputs or creating output. Use a new output directory.\n",
  "Reference/strand confirmation is deliberately unavailable here; evidence-supported declarations use the advanced native CLI.\n",
  sep = "")

if (sys.nframe() == 0L) {
  argv <- commandArgs(TRUE)
  if (identical(argv, "--help")) configured_analysis_help() else {
    dry <- "--dry-run" %in% argv
    if (sum(argv == "--dry-run") > 1L) stop("Duplicate --dry-run option.")
    argv <- argv[argv != "--dry-run"]
    if (!length(argv) || length(argv) %% 2L) stop("Use --config FILE [--action ACTION] [--dry-run], or --help.")
    flags <- argv[seq(1L, length(argv), 2L)]
    if (anyDuplicated(flags) || any(!flags %in% c("--config", "--action")) || !"--config" %in% flags) stop("Unknown, duplicate or missing configuration option.")
    args <- as.list(argv[seq(2L, length(argv), 2L)]); names(args) <- sub("^--", "", flags)
    args$dry_run <- dry
    do.call(run_configured_analysis, args)
  }
}
