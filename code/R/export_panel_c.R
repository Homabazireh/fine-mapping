#!/usr/bin/env Rscript
# Re-render only panel C from saved independent-R aggregate tables.
entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(entry)), "compare_haplotypes.R"))
argv <- commandArgs(TRUE)
if (length(argv) != 4L || !setequal(argv[c(1,3)], c("--run-dir", "--out"))) stop("Use: Rscript code/R/export_panel_c.R --run-dir COMPLETED_R_RUN --out NEW_DIRECTORY")
opts <- as.list(argv[c(2,4)]); names(opts) <- argv[c(1,3)]
run <- opts[["--run-dir"]]; out <- opts[["--out"]]
if (!file.exists(file.path(run, "COMPLETED.txt"))) stop("Incomplete R run.")
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) stop("Use a new output directory.")
summary <- jsonlite::fromJSON(file.path(run, "run_summary.json"), simplifyVector = FALSE)
args <- summary$parameters; args$out <- out; args$synthetic <- isTRUE(summary$synthetic)
cohort <- list(samples = seq_len(sum(unlist(summary$group_counts))), group_counts = unlist(summary$group_counts))
leads <- read_tsv(file.path(run, "lead_loci.tsv"))
loci <- lapply(seq_len(nrow(leads)), function(i) list(ident = leads$locus_id[i], chrom = leads$chrom[i], pos = as.numeric(leads$lead_pos[i])))
x <- read_tsv(file.path(run, "panel_c_data.tsv"))
for (k in c("lead_pos", "n_called", "n_group", "completeness")) x[[k]] <- suppressWarnings(as.numeric(x[[k]]))
data <- lapply(seq_len(nrow(x)), function(i) as.list(x[i, ]))
dir.create(out, recursive = TRUE)
draw_panel_c(args, cohort, loci, data)
file.copy(file.path(run, "panel_c_data.tsv"), file.path(out, "panel_c_data.tsv"))
cat("Saved only panel C from independent-R aggregate results: ", normalizePath(out), "\n", sep = "")
