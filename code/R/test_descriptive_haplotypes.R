#!/usr/bin/env Rscript
entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(entry)), "compare_haplotypes.R"))
source(file.path(dirname(normalizePath(entry)), "plot_descriptive_haplotypes.R"))
assert <- function(x, message) if (!isTRUE(x)) stop(message, call. = FALSE)
lead <- data.frame(locus_id = "L0001", chrom = "chr1", lead_pos = "100", status = "analysed")
summary <- list(group_counts = list(T2D = 10L, Prediabetes = 12L, Control = 14L), selected_loci = 1L, synthetic = TRUE, inputs_stable_during_read = TRUE)
rows <- list()
for (i in 1:3) for (h in 1:2) {
  n <- c(8, 10, 12)[i]; copies <- c(12, 10, 6)[i]; if (h == 2) copies <- 2*n-copies
  rows[[length(rows)+1L]] <- list(locus_id = "L0001", haplotype_id = paste0("L0001_H00", h), group = c("T2D", "Prediabetes", "Control")[i],
    n_group = c(10, 12, 14)[i], n_resolved = n, haplotype_copies = copies, resolved_chromosomes = 2*n, haplotype_frequency = copies/(2*n))
}
frequencies <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
data <- descriptive_haplotype_data(lead, frequencies, summary)
assert(length(data$plot) == 6L && length(data$differences) == 4L, "Haplotype patterns or requested contrasts omitted")
assert(data$differences[[1]]$descriptive_copy_frequency_difference == .5 && data$differences[[2]]$descriptive_copy_frequency_difference == .25, "Descriptive contrast direction wrong")
assert(data$interpretation[[1]]$haplotype_id_largest_absolute_descriptive_difference == "L0001_H001", "Largest-absolute tie did not use lexical haplotype ID")
assert(!any(c("p_value", "q_value", "CI95_low", "CI95_high") %in% names(data$differences[[1]])), "Inference attached to descriptive output")
unresolved <- frequencies
unresolved$n_resolved[unresolved$group == "Prediabetes"] <- 0
unresolved$haplotype_copies[unresolved$group == "Prediabetes"] <- 0
unresolved$resolved_chromosomes[unresolved$group == "Prediabetes"] <- 0
unresolved$haplotype_frequency[unresolved$group == "Prediabetes"] <- NA_real_
missing <- descriptive_haplotype_data(lead, unresolved, summary)
assert(all(vapply(Filter(function(r) r$contrast == "Prediabetes vs Control", missing$differences), function(r) is.na(r$descriptive_copy_frequency_difference), logical(1))), "Unresolved frequency became zero")
bad <- frequencies; bad$haplotype_copies[1] <- bad$haplotype_copies[1]+1
assert(inherits(try(descriptive_haplotype_data(lead, bad, summary), silent = TRUE), "try-error"), "Inconsistent copies accepted")
bad <- frequencies[-1, ]
assert(inherits(try(descriptive_haplotype_data(lead, bad, summary), silent = TRUE), "try-error"), "Missing haplotype/group row accepted")
bad <- rbind(frequencies, frequencies[1, ])
assert(inherits(try(descriptive_haplotype_data(lead, bad, summary), silent = TRUE), "try-error"), "Duplicate haplotype/group row accepted")
assert(inherits(try(descriptive_haplotype_data(rbind(lead, lead), frequencies, summary), silent = TRUE), "try-error"), "Duplicate lead IDs accepted")
extra <- summary; extra$group_counts$Unknown <- 5L
assert(inherits(try(descriptive_haplotype_data(lead, frequencies, extra), silent = TRUE), "try-error"), "Extra source group silently ignored")
assert(inherits(try(descriptive_haplotype_data(lead, frequencies[0, ], summary), silent = TRUE), "try-error"), "Empty data reported as available")
folder <- tempfile("native_R_descriptive_test_"); dir.create(folder)
draw_descriptive_haplotypes(data, summary, file.path(folder, "synthetic.pdf"), "pdf")
draw_descriptive_haplotypes(missing, summary, file.path(folder, "missing.pdf"), "pdf")
assert(all(file.info(list.files(folder, full.names = TRUE))$size > 1000), "Descriptive plot rendering failed")
unlink(folder, recursive = TRUE)
cat("PASS: native R descriptive denominator/copy conservation, all patterns, contrast signs, missingness, ties and no inference.\n")
