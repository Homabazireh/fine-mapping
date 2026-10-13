#!/usr/bin/env Rscript
# Verify that parallel execution preserves overlapping windows and missingness.
entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(entry)), "compare_haplotypes_parallel.R"))
assert <- function(x, message) if (!isTRUE(x)) stop(message)
loci <- list(list(chrom="chr1",start=80,end=120), list(chrom="chr1",start=110,end=150),
             list(chrom="chr1",start=480,end=520), list(chrom="chr2",start=80,end=120))
parts <- partition_windows(loci, 3L)
assert(identical(sort(unlist(parts)), seq_along(loci)), "Partition lost or duplicated a region")
assert(any(vapply(parts, function(x) all(c(1L,2L) %in% x), logical(1))), "Overlapping windows split across workers")
if (.Platform$OS.type == "windows" || !nzchar(Sys.which("bcftools"))) stop("Integration test requires Unix and bcftools")
folder <- tempfile("R_parallel_fixture_"); dir.create(folder)
samples <- paste0("synthetic_",1:6)
write.table(data.frame(sample=samples,group=rep(c("T2D","Prediabetes","Control"),each=2)),
  file.path(folder,"samples.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
writeLines(c("Chromsome\tPosition\tEffectAllele\tNonEffectAllele\tBeta\tSE\tPval\tliftover_strand",
  "1\t100\tG\tA\t0.2\t0.03\t1e-12\t+", "1\t130\tG\tA\t0.2\t0.03\t1e-11\t+",
  "1\t500\tG\tA\t0.2\t0.03\t1e-10\t+", "2\t25000000\tG\tA\t0.2\t0.03\t1e-9\t+"),file.path(folder,"gwas.tsv"))
header <- c("##fileformat=VCFv4.2", "##reference=GRCh38", "##contig=<ID=chr1,length=1000>", "##contig=<ID=chr2,length=100000000>",
 '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
 '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
 '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Quality">',
 '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
 paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",samples),collapse="\t"))
positions <- c(80,100,115,130,145,480,500,510,24999995,25000000,25000010)
chroms <- c(rep("chr1",8),rep("chr2",3))
records <- vapply(seq_along(positions),function(i) {
 calls <- c("0|0:20:40:1","0|1:20:40:1","1|1:20:40:1","0|1:20:40:1","0|0:20:40:1","1|1:20:40:1")
 if (positions[i]==115) calls[2] <- "./.:.:.:."
 paste(c(chroms[i],formatC(positions[i],format="f",digits=0),".","A","G",60,"PASS",".","GT:DP:GQ:PS",calls),collapse="\t")
},character(1))
vcf <- file.path(folder,"fixture.vcf"); compressed <- paste0(vcf,".gz")
writeLines(c(header,records),vcf)
assert(system2("bcftools",c("view","-Oz","-o",shQuote(compressed),shQuote(vcf)))==0,"Fixture compression failed")
assert(system2("bcftools",c("index",shQuote(compressed)))==0,"Fixture indexing failed")
args <- defaults(); args$vcf <- compressed; args$gwas <- file.path(folder,"gwas.tsv"); args$samples <- file.path(folder,"samples.tsv")
args$genome_build <- "GRCh38"; args$synthetic <- TRUE; args$figure <- FALSE; args$test_associations <- TRUE
args$top_loci <- 0; args$lead_gap <- 10; args$flank <- 20
args$out <- file.path(folder,"sequential"); args$workers <- 1L; sequential <- run_analysis(args)
args$out <- file.path(folder,"parallel"); args$workers <- 3L; parallel <- run_analysis(args)
for (name in union(names(sequential$scan_counts),names(parallel$scan_counts)))
 assert(identical(as.numeric(sequential$scan_counts[[name]]),as.numeric(parallel$scan_counts[[name]])),paste("Counter mismatch",name))
assert(sequential$scan_counts$records_in_windows==length(records),"Overlapping windows counted duplicate records")
keys <- read_tsv(file.path(args$out,"regional_variant_frequencies.tsv.gz"))$variant_id
assert("chr2:25000000:A:G" %in% keys && !any(grepl("e\\+", keys)), "Genomic identifiers used scientific notation")
for (name in c("lead_loci.tsv","regional_variant_frequencies.tsv.gz","regional_variant_qc.tsv.gz","selected_markers.tsv",
 "sample_haplotypes.tsv","haplotype_group_frequencies.tsv","haplotype_comparisons.tsv","haplotype_overall_comparisons.tsv",
 "group_phase_qc.tsv","selected_marker_genotypes.tsv.gz","panel_c_data.tsv")) {
 a <- read_tsv(file.path(folder,"sequential",name)); b <- read_tsv(file.path(folder,"parallel",name))
 if(nrow(a)) a <- a[do.call(order,a),,drop=FALSE]
 if(nrow(b)) b <- b[do.call(order,b),,drop=FALSE]
 rownames(a)<-rownames(b)<-NULL
 assert(identical(a,b),paste("Scientific table differs",name))
}
# A failed worker must not produce a completed run or discard diagnostics.
original_scanner <- sequential_scan_vcf
sequential_scan_vcf <- function(...) stop("Injected worker extraction failure")
args$out <- file.path(folder,"worker_failure")
failure <- suppressWarnings(tryCatch(run_analysis(args), error=identity))
sequential_scan_vcf <- original_scanner
assert(inherits(failure,"error"),"Worker failure was silently accepted")
assert(!file.exists(file.path(args$out,"COMPLETED.txt")),"Failed worker yielded a completed analysis")
assert(dir.exists(file.path(args$out,"scan_workers")),"Worker failure discarded diagnostics")
unlink(folder,recursive=TRUE)
cat("PASS: native parallel R matches sequential R across 11 tables and all scan counters; overlapping windows, chromosome boundaries and missing calls preserved.\n")
