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

## Figure generation

```bash
.venv/bin/python code/export_panel_c.py \
  --run-dir results/all_regions_run --out results/all_regions_run/panel_c
.venv/bin/python code/plot_haplotype_contrasts.py \
  --run-dir results/all_regions_run --out results/all_regions_run/contrasts
.venv/bin/python code/plot_all_regions.py \
  --run-dir results/all_regions_run --out results/all_regions_run/regions --workers 4
.venv/bin/python code/plot_descriptive_haplotypes.py \
  --run-dir results/all_regions_run --out results/all_regions_run/descriptive_haplotypes
```

Panel C displays lead-genotype completeness. It must not be described as a haplotype frequency plot. The contrast figure displays eligible carrier-frequency differences with confidence intervals and adjusted p-values, or an explicit reason why an estimate is unavailable.

The all-region plotting script draws a separate picture for each region using every available passing biallelic SNP in its saved window. It records omitted indels and multiallelic sites explicitly and preserves the full regional allele table. Completeness is called participants divided by the full group size. Regional ALT frequency is ALT copies divided by called chromosomes; it is conditional on successful calls and can be biased when calls are selective. Missing frequencies remain unavailable. Regional scatterplots do not establish phase, shared ancestry, association or causality.

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
