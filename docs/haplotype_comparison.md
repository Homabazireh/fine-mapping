# Reproduce the three-group comparison

The final participant table defines 49 Diabetes, 48 Prediabetes and 50 Control
participants. Use `data/raw/participant_groups_final.tsv` and the verified alias
map; earlier participant lists are superseded. Inputs use GRCh38.

## Python workflow

Run from the repository root with the original data restored locally:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-analysis.txt
.venv/bin/python code/make_sample_template.py \
  --metadata data/raw/participant_groups_final.tsv \
  --id-map data/intermediate/participant_review/verified_id_map.tsv \
  --out data/intermediate/samples_reproduced.tsv
.venv/bin/python code/compare_haplotypes.py \
  --samples data/intermediate/samples_reproduced.tsv \
  --out results/reproduced/cohort --genome-build GRCh38 \
  --test-associations --no-plot
.venv/bin/python code/audit_finemapping.py \
  --run-dir results/reproduced/cohort \
  --gwas data/raw/EUR_diagram_t2d_1e-8.tsv \
  --legacy-dir ../GWAS_regions/SuSiE_all_leads \
  --out results/reproduced/variant_evidence
.venv/bin/python code/publish_haplotype_comparison.py \
  --run-dir results/reproduced/cohort \
  --evidence-dir results/reproduced/variant_evidence \
  --out results/reproduced/publication
```

Use new output paths for each run. Only panel C is rendered, in PNG/PDF/SVG.
The private data release restores the original `Fine_mapping/new_analysis`
layout so that the example relative paths resolve. A public checkout includes
the code and aggregate publication, but requires local data for a cohort rerun.

The completed 13 October run is under
`results/haplotype_comparison_20261013/`. The publication is copied to
`docs/figures/`; its tables have no participant identifiers or individual
genotypes. Local detailed run files must not be staged into the public repo.

## Statistical interpretation

- Chromosome-copy frequencies are descriptive and use twice the number of
  resolved participants as their denominator. Carrier frequencies count each
  participant once, including participants with two copies of the haplotype.
- Each haplotype receives an exact three-group Fisher–Freeman–Halton carrier
  test and pairwise two-sided Fisher carrier tests when eligible. Pairwise
  effects are conditional odds ratios with exact 95% intervals and carrier
  frequency differences, group A minus group B, with Newcombe–Wilson intervals.
- BH correction is applied across all eligible locus/haplotype tests separately
  for the overall and pairwise families. Intervals are pointwise 95% intervals,
  not multiplicity-adjusted intervals. Positive differences mean higher carrier
  frequency in group A. Tests are exploratory and unadjusted for covariates.
- Inference additionally requires verified reference-call provenance, at least
  10 resolved participants and 80% resolution in every compared group, and at
  least three pooled carriers and noncarriers. `--test-associations` requests
  tests but does not attest genotype provenance. Only use
  `--reference-genotypes-confirmed` after verifying genuine reference calls.
- No phased haplotypes are stitched across incompatible phase sets. Missing
  calls remain missing. All groups use the same phenotype-blind marker panel.

Current results contain zero eligible panels and zero tests. Header-only
frequency/association/candidate tables mean unavailable, not equal frequencies,
zero effect, a p-value of one, or an absence of causal variants.

## Candidate variants and existing fine-mapping

The evidence audit joins the 20 lead variants to the original external GWAS
with exact allele matching and preserves extremely small p-values as text.
Those are external association signals, not associated cohort haplotypes.

Legacy SuSiE runs are audited for genotype availability, allele/z-score
alignment, sample-size information, model status and recorded diagnostics.
Their exploratory PIPs are not republished as validated causal rankings.
Convergence alone does not establish appropriate LD or valid genotype calls.
The audit reports its source hashes and specific reasons in aggregate tables.

To resolve candidate variants, provide reliable joint genotypes including
reference calls and phase; validate group missingness, ancestry, relatedness
and covariates; then use harmonized full regional GWAS statistics and suitable
ancestry-matched LD with fine-mapping diagnostics. Independent replication and
functional evidence would strengthen causal interpretation.

Methods: [R Fisher-test documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/fisher.test.html),
[Newcombe–Wilson implementation reference](https://www.statsmodels.org/stable/_modules/statsmodels/stats/proportion.html),
[SuSiE-RSS diagnostics](https://stephenslab.github.io/susieR/articles/susierss_diagnostic.html),
[SuSiE-RSS API](https://stephenslab.github.io/susieR/reference/susie_rss.html).

## R and validation

The [R README](../code/R/README.md) documents the independent native R workflow,
dependencies, real-cohort run, synthetic regression checks and parity checks.
R does not call the Python implementation to calculate its analysis results.

```bash
.venv/bin/python code/run_tests.py
git diff --check
```

Synthetic fixtures test valid and invalid phase, genotype quality, reference
provenance, empty outputs, exact tests, uncertainty intervals and correction
families. Synthetic results are never used in the cohort publication.
