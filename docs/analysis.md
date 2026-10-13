# Analysis workflow

The public GWAS selects regions; the VCF supplies genotypes; verified participant metadata supplies Diabetes, Prediabetes and Control labels. Results and figures must remain private. Run commands from the repository root and use a new output directory for every analysis.

## Inputs and preparation

Configure paths in `config/settings.json`. Required inputs are a VCF, GWAS table and a `sample`/`group` TSV with exact VCF IDs. Use the final authoritative metadata table and independently verified one-to-one aliases; never infer an identity match from a similar identifier. The preparation command is in the repository README. Check that all inputs already use the declared genome build: the workflow does not perform liftover.

The GWAS reader preserves extremely small p-values using decimal ranking. It selects variants at the configured significance threshold, groups consecutive positions by the configured gap and uses the lowest-p variant as each region's lead. This is physical-distance selection, not LD clumping. The default subset is controlled by `top_loci`; `--top-loci 0` selects every qualifying region.

## Full-region analysis

```bash
.venv/bin/python code/compare_haplotypes.py \
  --samples data/intermediate/samples.tsv --genome-build GRCh38 \
  --top-loci 0 --test-associations --out results/all_regions_run
```

The VCF is scanned once. Regional tables retain the qualifying SNP/indel records in each window; a lead is an anchor and does not replace the surrounding markers. Overlapping windows retain each region's own membership. Genotype DP/GQ, PASS, pooled completeness, allele orientation and phase checks remain active. Never convert an absent call to a homozygous-reference genotype or reduce thresholds to force a comparison.

Marker panels use the same variants in every group. A participant with more than one heterozygous marker requires consistent, explicit phase-set information. Descriptive haplotype frequencies use resolved chromosome copies. Inferential comparisons use participant carriers, with reference-genotype provenance and minimum resolution/count checks. See [the statistical specification](haplotype_comparison.md).

To diagnose absent reference calls independently of genotype-quality filtering, audit raw VCF genotypes across the union of selected windows:

```bash
.venv/bin/python code/audit_reference_genotypes.py \
  --vcf data/raw/full_LR-WGS.vcf.gz \
  --leads results/all_regions_run/lead_loci.tsv \
  --out results/reference_genotype_audit_new
```

This bounded-memory diagnostic needs `bcftools` and `awk`. It records aggregate raw genotype categories, missing DP/GQ at matching leads, source checksums and merge provenance, without participant identifiers. Its results remain private. It diagnoses unavailable reference information; it cannot reconstruct genotypes or authorize the reference-provenance flag.

## Figure generation

```bash
.venv/bin/python code/export_panel_c.py \
  --run-dir results/all_regions_run --out results/all_regions_run/panel_c
.venv/bin/python code/plot_haplotype_contrasts.py \
  --run-dir results/all_regions_run --out results/all_regions_run/contrasts
.venv/bin/python code/plot_all_regions.py \
  --run-dir results/all_regions_run --out results/all_regions_run/regions \
  --include-indels --workers 4
.venv/bin/python code/plot_descriptive_haplotypes.py \
  --run-dir results/all_regions_run --out results/all_regions_run/descriptive_haplotypes
```

Panel C displays lead-genotype completeness. It must not be described as a haplotype frequency plot. The contrast figure displays observed participant-carrier differences for every resolved haplotype, even when inference is withheld. Its overview contains all available panels, with separate regional pictures and interpretations linked from `index.md`; the complete table retains excluded regions and their reasons. Descriptive observations have no confidence intervals or significance marks. Valid tested comparisons additionally display their confidence intervals and original full-family adjusted p-values in the regional pictures. Missing observations remain unavailable, while a genuine observed zero is shown at zero.

Interpretation reports count identified-haplotype comparisons separately from region-status placeholders. An excluded region contributes two requested Control-contrast placeholders without an identified haplotype; these are not two failed haplotype tests. Regional reports state their scope within the full source run. The first reported inference gate may hide additional participant-resolution or carrier-count failures, and a post-QC reference-call count cannot establish what was present before QC. Use the independent raw-genotype audit for that distinction. Contrast indexes and manifests record the complete source directory.

The all-region plotting script draws a separate picture for each region. With `--include-indels`, it includes passing biallelic SNPs and explicit small-indel sequences within the saved analysis's allele-length limit. Circles identify SNPs; triangles identify indels; colors identify participant groups. Tables and manifests include `variant_type`, separate SNP/indel counts and `plotted_variant_count`. Multiallelic and unsupported records have explicit omission counts and remain in the full regional allele table. Omitting the option preserves the earlier SNP-only schema.

Completeness is called participants divided by the full group size. Regional ALT frequency is ALT copies divided by called chromosomes; it is conditional on successful calls and can be biased when calls are selective. A zero-call site has zero completeness and unavailable ALT frequency. Allele-copy counts are counts of the two diploid alleles, not copy-number-variation estimates. Regional scatterplots do not establish phase, shared ancestry, association or causality.

When descriptive haplotype panels exist but inference remains unavailable, `plot_descriptive_haplotypes.py` shows every observed local pattern: three chromosome-copy frequencies and two descriptive disease-minus-Control differences. It recalculates fractions from counts and reports resolution denominators. The picture does not assign p-values, confidence intervals or significance to these untested differences.

The optional report workflow is `code/build_report.py --run-dir RUN --no-preflight --output-dir OUTPUT --render`; rendering needs Quarto. Reports and embedded images are private outputs. The original example-region selector in the legacy report is a layout example; use the dedicated all-region workflow for exhaustive regional exports.

## Outputs and evidence

| File | Meaning |
| --- | --- |
| `lead_loci.tsv` | Region boundaries, leads and eligibility or exclusion reasons |
| `regional_variant_frequencies.tsv.gz` | Regional allele counts, called denominators and group totals |
| `regional_variant_qc.tsv.gz` | Variant-level quality and marker eligibility |
| `selected_markers.tsv` | Common marker panels |
| `sample_haplotypes.tsv` | Participant-level haplotypes or unresolved status; private data |
| `haplotype_group_frequencies.tsv` | Descriptive chromosome-copy frequencies and carrier counts |
| `group_phase_qc.tsv` | Resolved participants and exclusions by group |
| `haplotype_overall_comparisons.tsv` | Overall carrier tests and adjusted p-values |
| `haplotype_comparisons.tsv` | Pairwise effects, intervals and adjusted p-values |
| `run_summary.json` | Input provenance, parameters, counts and warnings |

A completed run records `COMPLETED.txt`. Keep its source tables and provenance together. Genomic TSV coordinates are one-based; BED start coordinates are zero-based and the end is exclusive. The independent [R workflow](../code/R/README.md) can reproduce analysis tables and figures, with `verify_parity.R` checking numerical agreement.

`code/audit_finemapping.py` and its R equivalent audit existing fine-mapping evidence, including harmonization and LD limitations. Valid candidate localization requires suitable genotypes, full harmonized regional summary statistics, sample-size information, appropriate LD and fine-mapping diagnostics. Functional evidence and independent replication strengthen causal interpretation. Unavailable estimates do not demonstrate equal group frequencies.

For a separate exploratory model using documented public GWAS statistics and independent reference genotypes, follow the [external SNP/indel extraction, LD and SuSiE-RSS workflow](external_finemapping.md). Its outputs describe the external study and remain separate from cohort haplotype tests.
