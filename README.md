# Fine mapping — reproducible analysis code

Compare local haplotypes across Diabetes, Prediabetes and Control using a shared, phenotype-blind marker panel. Python and independent R workflows prepare data, apply genotype and phase checks, calculate eligible group comparisons, and draw the requested figures.

This public repository contains code, tests and methodological documentation only. All cohort data, results, figures and interpretations belong in the [private fine-mapping-2 repository](https://github.com/Homabazireh/fine-mapping-2). Access to that repository is required to view or restore them. Removing outputs here does not remove their earlier public Git history.

## Reproduce an analysis

Use Python 3.10+ and R 4.3+ with `jsonlite`; install `bcftools` for the R VCF reader. Supply the inputs documented in [data/raw/README.md](data/raw/README.md), or restore the private release according to its `PRIVATE_DATA.md`.

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-analysis.txt
.venv/bin/python code/make_sample_template.py \
  --metadata data/raw/participant_groups_final.tsv \
  --id-map data/intermediate/participant_review/verified_id_map.tsv \
  --out data/intermediate/samples_for_review.tsv
.venv/bin/python code/compare_haplotypes.py \
  --samples data/intermediate/samples_for_review.tsv \
  --genome-build GRCh38 --test-associations \
  --out results/cohort_run
.venv/bin/python code/export_panel_c.py \
  --run-dir results/cohort_run --out results/cohort_run/panel_c
.venv/bin/python code/plot_haplotype_contrasts.py \
  --run-dir results/cohort_run --out results/cohort_run/contrasts
```

Only declare `--reference-genotypes-confirmed` or `--gwas-alleles-forward-confirmed` after verifying the corresponding evidence. A request for association tests does not bypass eligibility checks. Use `--top-loci 0` to analyze every GWAS-selected region rather than the default subset; preserve a new output folder for each run.

[Analysis instructions](docs/analysis.md) describe inputs, all-region plots and reports. [Statistical methods](docs/haplotype_comparison.md) explain effect sizes and uncertainty. [R instructions](code/R/README.md) provide the independent implementation and parity checks. [The revised task](docs/request.md) records the requested scope.

## Validation

```bash
.venv/bin/python code/run_tests.py
Rscript code/R/test_analysis.R
Rscript code/R/test_audit.R
git diff --check
```

Tests use synthetic fixtures and do not establish results for a real cohort. The analysis preserves missing calls, validates phase, distinguishes genotype completeness from haplotype frequency, and reports unavailable estimates explicitly. Association and fine-mapping evidence do not by themselves establish causality.

See [AGENTS.md](AGENTS.md) for scientific and publishing requirements and [GitHub instructions](docs/github.md) for the separation between public code and private outputs.
