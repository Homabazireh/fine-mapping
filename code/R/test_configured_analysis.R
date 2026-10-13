#!/usr/bin/env Rscript
# Synthetic inputs only; exercise the configured entry point and child processes.
entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(entry)), "run_configured_analysis.R"))
assert <- function(value, message) if (!isTRUE(value)) stop(message)
fails <- function(expr, pattern) {
  error <- tryCatch({force(expr); NULL}, error = identity)
  assert(inherits(error, "error") && grepl(pattern, conditionMessage(error)), paste("Expected error:", pattern))
}
root <- dirname(dirname(R_WORKFLOW_DIRECTORY))
folder <- tempfile("R_workflow_"); dir.create(folder)
config_folder <- file.path(folder, "config with spaces"); dir.create(config_folder)
config <- file.path(config_folder, "workflow.json")
save_config <- function(x) jsonlite::write_json(x, config, auto_unbox = TRUE, pretty = TRUE, null = "null")
value_after <- function(command, flag) command$args[match(flag, command$args) + 1L]

cfg <- list(action = "regional_plots", actions = list(
  regional_plots = list(input = list(run_dir = "../missing source"), out = "../new plots", options = list(include_indels = TRUE)),
  cohort_analysis = list(input = list(config = file.path(root, "config/settings.json")), out = "../cohort", options = list(genome_build = "GRCh38", workers = 1)),
  descriptive_haplotypes = list(input = list(run_dir = "../missing source"), out = "../descriptive", options = list(columns = 2)),
  haplotype_contrasts = list(input = list(run_dir = "../missing source"), out = "../contrasts", options = list()),
  external_finemapping = list(input = list(config = "external.json"), out = "../external new", options = list(max_effects = 4, coverage = .9, workers = 2, variant_types = "both"))))
save_config(cfg)
# A sourced file must not interpret the caller's CLI or start a workflow.
assert(is.function(run_configured_analysis) && !dir.exists(file.path(folder, "new plots")), "source() started an analysis")
old <- getwd(); setwd(tempdir())
command <- run_configured_analysis(config, dry_run = TRUE)
setwd(old)
assert(identical(basename(command$script), "plot_all_regions.R") && "--include-indels" %in% command$args, "SNP+indel plotting not selected")
assert(identical(value_after(command, "--run-dir"), file.path(config_folder, "../missing source")), "Paths resolved against working directory instead of config")
assert(!dir.exists(file.path(folder, "new plots")), "Dry run created outputs")
cfg$actions$descriptive_haplotypes$options$columns <- "invalid while inactive"
save_config(cfg)
invisible(configured_analysis_command(config)) # only selected settings are validated
fails(configured_analysis_command(config, "descriptive_haplotypes"), "finite number")
cfg$actions$descriptive_haplotypes$options$columns <- 2
cfg$actions$regional_plots$options$include_indels <- "true"
save_config(cfg); fails(configured_analysis_command(config), "true or false")
cfg$actions$regional_plots$options$include_indels <- FALSE
save_config(cfg); assert(!"--include-indels" %in% configured_analysis_command(config)$args, "False flag became enabled")
cfg$actions$regional_plots$options$include_indels <- TRUE
cfg$actions$regional_plots$options$executable <- "untrusted"
save_config(cfg); fails(configured_analysis_command(config), "Unknown options")
cfg$actions$regional_plots$options$executable <- NULL
save_config(cfg)
command <- configured_analysis_command(config, "cohort_analysis")
assert(basename(command$script) == "compare_haplotypes.R", "Sequential cohort action uses wrong script")
assert(all(c("--no-reference-genotypes-confirmed", "--no-gwas-alleles-forward-confirmed") %in% command$args), "Workflow enabled scientific provenance declarations")
cfg$actions$cohort_analysis$options$workers <- 2
save_config(cfg)
assert(basename(configured_analysis_command(config, "cohort_analysis")$script) == "compare_haplotypes_parallel.R", "Parallel cohort action not selected")
cfg$actions$cohort_analysis$options$reference_genotypes_confirmed <- TRUE
save_config(cfg); fails(configured_analysis_command(config, "cohort_analysis"), "Unknown options")
cfg$actions$cohort_analysis$options$reference_genotypes_confirmed <- NULL
cfg$actions$cohort_analysis$options$workers <- 1.5
save_config(cfg); fails(configured_analysis_command(config, "cohort_analysis"), "integer")
cfg$actions$cohort_analysis$options$workers <- 1
save_config(cfg)
fails(configured_analysis_command(config, "external_finemapping"), "does not exist")
jsonlite::write_json(list(paths = list(input_dir = "../prepared", out = "../external"), analysis = list(max_effects = 10)), file.path(config_folder, "external.json"), auto_unbox = TRUE)
command <- configured_analysis_command(config, "external_finemapping")
assert(value_after(command, "--config") == file.path(config_folder, "external.json") && value_after(command, "--max-effects") == "4" && value_after(command, "--coverage") == "0.9", "External config or options were lost")

# Verify literal argv round-tripping through the actual system2 shell boundary.
echo_script <- file.path(folder, "echo arguments.R")
capture_file <- file.path(folder, "captured.rds")
writeLines("args <- commandArgs(TRUE); saveRDS(args[-1], args[1])", echo_script)
literal <- c("spaces and 'single quotes' and \"double quotes\"", "$(touch SHOULD_NOT_EXIST)", "`touch ALSO_NOT_CREATED`", "semi;colon", "a&b", "percent%sign")
probe <- list(executable = file.path(R.home("bin"), "Rscript"), script = echo_script, args = c(capture_file, literal))
setwd(folder); workflow_invoke(probe); setwd(old)
assert(identical(readRDS(capture_file), literal), "Arguments were interpreted by the shell")
assert(!file.exists(file.path(folder, "SHOULD_NOT_EXIST")) && !file.exists(file.path(folder, "ALSO_NOT_CREATED")), "Shell substitution executed")
writeLines("quit(status=7)", echo_script)
fails(workflow_invoke(probe), "exit status 7")

# End-to-end: an unindexed synthetic VCF, cohort analysis, and all three plots.
inputs <- file.path(folder, "synthetic inputs 'quoted'"); dir.create(inputs)
samples <- paste0("synthetic_", 1:6)
write.table(data.frame(sample = samples, group = rep(c("T2D", "Prediabetes", "Control"), each = 2)),
  file.path(inputs, "samples.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c("Chromsome\tPosition\tEffectAllele\tNonEffectAllele\tBeta\tSE\tPval\tliftover_strand",
  "1\t100\tG\tA\t0.2\t0.03\t1e-12\t+"), file.path(inputs, "gwas.tsv"))
header <- c("##fileformat=VCFv4.2", "##reference=GRCh38", "##contig=<ID=chr1,length=1000>",
  paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", samples), collapse = "\t"))
gt <- c("0|0:20:40:1", "0|1:20:40:1", "1|1:20:40:1", "0|1:20:40:1", "0|0:20:40:1", "1|1:20:40:1")
records <- vapply(1:2, function(i) paste(c("chr1", c(100, 110)[i], ".", c("A", "AT")[i], c("G", "A")[i], "60", "PASS", ".", "GT:DP:GQ:PS", gt), collapse = "\t"), character(1))
writeLines(c(header, records), file.path(inputs, "cohort.vcf"))
cfg$actions$cohort_analysis$input <- list(config = file.path(root, "config/settings.json"), vcf = file.path(inputs, "cohort.vcf"),
  gwas = file.path(inputs, "gwas.tsv"), samples = file.path(inputs, "samples.tsv"))
cfg$actions$cohort_analysis$options <- list(genome_build = "GRCh38", workers = 1, top_loci = 0, synthetic = TRUE, figure = FALSE, test_associations = TRUE)
for (action in c("regional_plots", "haplotype_contrasts", "descriptive_haplotypes")) cfg$actions[[action]]$input$run_dir <- "../cohort"
save_config(cfg)
run_configured_analysis(config, "cohort_analysis")
summary <- jsonlite::fromJSON(file.path(folder, "cohort/run_summary.json"))
assert(summary$synthetic && summary$haplotype_carrier_tests_performed == 0 && !summary$parameters$reference_genotypes_confirmed, "Wrapper weakened missing-provenance inference gate")
for (action in c("regional_plots", "haplotype_contrasts", "descriptive_haplotypes")) run_configured_analysis(config, action)
assert(file.exists(file.path(folder, "new plots/index.md")) && file.exists(file.path(folder, "contrasts/INTERPRETATION.txt")) && file.exists(file.path(folder, "descriptive/descriptive_haplotypes.png")), "A configured plotting action produced no index/figure")
assert(file.exists(file.path(folder, "new plots/L0001_chr1_100/region.png")), "Configured SNP+indel regional picture absent")
fails(run_configured_analysis(config, "cohort_analysis"), "exit status")
assert(file.exists(file.path(folder, "cohort/COMPLETED.txt")), "Repeat execution damaged the completed analysis")
unlink(folder, recursive = TRUE)
cat("PASS: configured R workflow, relative paths, selected actions, provenance gates, literal shell arguments, error propagation and synthetic cohort/plot integration.\n")
