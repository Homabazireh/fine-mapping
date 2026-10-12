# ZFHX3 haplotype comparison and its relevance to our diabetes cohort

Reviewed 12 October 2026 against the completed `cohort_no_preflight` run.

**We have not performed the same analysis as this paper.** We share the idea of comparing nearby DNA variants on individual chromosome copies. The paper investigates a shared disease-associated chromosome background across families. Our intended analysis compares how frequently local patterns occur in Diabetes, Prediabetes, and Control participants. Our completed run produced no usable haplotype panels, so those group comparisons are still unavailable.

## The paper

Figueroa, K. P., Gross, C., Buena-Atienza, E., et al. (2024). *A GGC-repeat expansion in ZFHX3 encoding polyglycine causes spinocerebellar ataxia type 4 and impairs autophagy.* **Nature Genetics 56**, 1080–1089. Published online 29 April 2024. DOI: [10.1038/s41588-024-01719-5](https://doi.org/10.1038/s41588-024-01719-5). [PubMed record](https://pubmed.ncbi.nlm.nih.gov/38684900/).

This is a study of **spinocerebellar ataxia type 4 (SCA4)**, a neurological disorder. The authors identified a GGC repeat expansion in **ZFHX3** using long-read sequencing. Rare nearby variants supported a distant shared founder origin. The paper also reports cellular evidence involving impaired autophagy. These findings do not establish a diabetes mechanism. [Article abstract](https://www.nature.com/articles/s41588-024-01719-5).

## What the haplotype comparison means

A **haplotype** is a pattern of DNA variants on one chromosome copy. Most people have two copies of each autosome. **Phasing** determines which variants belong together on each copy; having a list of variants alone does not always tell us this.

In **Figure 1d**, the authors compare single-base variants within **50 kb on either side** of the repeat expansion. The expansion-bearing chromosome backgrounds resemble one another in cases from Utah, Lübeck, and Tübingen. Only the Utah and Lübeck cases share all six highlighted rare variants. In plain language, the question is: **“Do these affected families carry a similar stretch of DNA around the expansion?”** This contributes to the founder interpretation; visual similarity by itself is not a complete causal demonstration. [Figure 1 and its caption](https://www.nature.com/articles/s41588-024-01719-5/figures/1).

## Comparison with our code and results

| Aspect | Paper and associated study code | Our cohort analysis |
| --- | --- | --- |
| Main question | Shared chromosome background around the SCA4 expansion | Frequencies of local patterns across Diabetes, Prediabetes, and Control |
| Region | ZFHX3 repeat region | 20 regions selected from public type 2 diabetes GWAS results |
| Local window | Figure 1d examines ±50 kb | ±50 kb around each selected GWAS lead |
| Variants | Nearby rare variants and the repeat expansion | Lead plus up to seven eligible nearby SNP/small-indel markers |
| Phasing | Public study workflow includes LongPhase from reads and variants | Reads existing VCF genotype and phase-set fields; performs no new read-based or statistical phasing |
| Repeat sizing | Study workflow includes TRGT | No repeat-expansion measurement |
| Current outcome | Published comparison of expansion-associated backgrounds | Zero usable panels and zero haplotype association tests in our completed run |

The study's public workflow includes [LongPhase phasing and haplotagging](https://github.com/imgag/sca4_project/blob/876195701730e4daca1cc5bdd380ccbaf2698380/rules/phasing.smk), [joint genotyping from gVCFs using GLnexus](https://github.com/imgag/sca4_project/blob/876195701730e4daca1cc5bdd380ccbaf2698380/rules/merge_variants.smk), and [TRGT repeat genotyping](https://github.com/imgag/sca4_project/blob/876195701730e4daca1cc5bdd380ccbaf2698380/rules/trgt.smk). These are study-code capabilities; they do not specify every step used for the published Figure 1d.

Our [comparison code](../../code/compare_haplotypes.py) requires consistent existing phase sets when a participant has two or more heterozygous markers. With zero or one heterozygous marker, the unordered chromosome pair can be resolved without additional phasing. It gives identical patterns the same local haplotype ID, using the same marker panel for every group. Such matching does not establish common ancestry. Region selection uses physical distance, not linkage-disequilibrium clumping or statistical fine-mapping.

Our code already includes a `plot_locus` function with **two chromosome rows per person and the same marker columns across groups**. However, plotting was disabled in this run, and no usable panels existed to plot. The cohort completeness heatmap in the earlier report is a different figure: its colors show the fraction of usable lead genotypes, not shared haplotypes or diabetes risk.

## What your data currently show

The completed direct-run summary (`results/cohort_no_preflight/run_summary.json` in the local analysis workspace) records **147 participants: 49 Diabetes, 48 Prediabetes, and 50 Control**. All participants have assigned groups; there are no unresolved group assignments. The [published panel c and its report](../../docs/figures/README.md) display the lead-genotype completeness measurements.

Of 20 regions, 19 were recorded as failing the 90% usable-lead-genotype requirement. One was recorded as having unconfirmed palindromic allele orientation; it also falls below the completeness requirement. Consequently, **there are no estimated haplotype frequencies or tested group differences**. This does not mean that the groups have identical haplotypes.

The run also found no high-quality homozygous reference (`0/0`) calls in the selected windows. This warrants checking how the multi-sample VCF was generated. It does not prove a particular processing error. Missing calls must remain missing unless genuine genotype evidence resolves them. Skipping the separate VCF precheck did not remove the analysis's internal requirements. The full local explanation and report are saved as `results/cohort_no_preflight/RESULTS.txt` and `results/cohort_no_preflight/report.html`.

## How the paper can guide the next analysis

1. **Establish usable genotypes and phase.** Check the source of missing and reference calls, and obtain reliable calls from the available per-sample data where needed. Determine which variants are on the same chromosome copy. Lowering a completeness threshold alone does not reconstruct missing genotypes or phase.
2. **Display actual chromosome patterns at selected loci.** Use the same marker positions for all three groups, mark the GWAS lead, and show unresolved patterns explicitly. The existing plot function can provide an initial view once panels are available. A denser, physically scaled display would better support examination of local sharing; eight markers cannot define an ancestral chromosome segment.
3. **Compare frequencies with clear denominators.** Report the number of resolved participants in each group and how often each pattern occurs. Any inferential analysis needs a suitable plan for missingness, ancestry, relatedness, relevant covariates, and multiple comparisons.

These are proposed extensions. No new phasing, repeat analysis, founder analysis, relaxed-filter run, or association test was performed while adding this literature review. A family-based founder study would require a separate design and supporting data.

## Saved sources and access

- [Publisher supplementary information PDF](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41588-024-01719-5/MediaObjects/41588_2024_1719_MOESM1_ESM.pdf): Supplementary Figures 1–5 and Tables 1–2. The downloaded PDF remains a local reading copy.
- [BibTeX citation](../references.bib).
- [Source URLs, repository revision, and download checksum](sources.json).

The main article PDF was unavailable through the public publisher access used here. The comparison above relies on the public abstract, Figure 1 caption, linked study code, and our local code/results; it is not a full-methods replication of the paper.
