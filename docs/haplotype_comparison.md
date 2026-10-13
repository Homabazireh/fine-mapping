# Haplotype comparison methods

Use a shared phenotype-blind marker panel within each region and preserve missing genotypes. Existing VCF phase sets determine whether a participant's two haplotypes can be resolved. Cohort findings and interpretation reports are stored only in the private repository.

## Frequency and inference

Haplotype chromosome-copy frequency is descriptive: copies of a pattern divided by all resolved chromosome copies in that group. The inferential unit is the participant, classified as carrier or noncarrier of a particular haplotype. Do not treat the two chromosomes from one participant as independent observations.

For eligible hypotheses, the overall three-group comparison uses an exact Fisher–Freeman–Halton test of the 3×2 carrier table. Pairwise comparisons use two-sided Fisher exact tests. Report conditional odds ratios with exact 95% confidence intervals and carrier-frequency differences (group A minus group B) with Newcombe–Wilson 95% confidence intervals. A positive difference means more carriers in group A. Intervals are pointwise, not multiplicity-adjusted.

Apply Benjamini–Hochberg correction separately across all tested overall hypotheses and across all tested pairwise hypotheses. The pairwise family includes all three group pairs, even if a figure highlights only the two control comparisons. Do not recompute correction over only displayed contrasts.

Inference requires verified reference-call provenance, at least ten resolved participants in each group, at least 80% resolution in each group, and at least three pooled carriers and noncarriers. A requested test can remain unavailable when eligibility fails. These are exploratory, unadjusted tests; covariate- or relatedness-aware inference needs a separately validated design.

## Figures and causal interpretation

Panel C is a genotype-completeness display. The haplotype contrast figure shows eligible carrier differences and uncertainty; unavailable effects must not be placed at zero. Regional marker pictures include all qualifying regional SNPs, report called denominators and distinguish conditional ALT frequency from cohort frequency. They cannot resolve haplotypes or prove a group difference by visual inspection.

An associated haplotype can help define candidate regions but does not establish a causal variant. Audit reference-call completeness, allele harmonization, summary-statistic coverage, sample size, LD ancestry and model diagnostics before using fine-mapping ranks. Replication, segregation where relevant, molecular evidence and appropriate functional experiments can further distinguish candidate mechanisms. Record unsupported conclusions as unresolved.
