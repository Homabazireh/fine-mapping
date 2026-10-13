# Fine mapping — reproducible analysis code

Compare local SNP/small-indel haplotypes across Diabetes, Prediabetes and Control using a shared, phenotype-blind marker panel. Python and independent R workflows prepare data, apply genotype and phase checks, calculate eligible group comparisons, and draw the requested figures. Allele-copy counts refer to diploid genotypes; this is not a copy-number-variation caller.

This public repository contains code, tests and methodological documentation only. All cohort data, results, figures and interpretations belong in the [private fine-mapping-2 repository](https://github.com/Homabazireh/fine-mapping-2). Access to that repository is required to view or restore them. Removing outputs here does not remove their earlier public Git history.

## Reproduce an analysis

Use Python 3.10+ and R 4.3+ with `jsonlite`; install `bcftools` for the R VCF reader. External LD preparation additionally uses NumPy, and external fine-mapping requires the R packages `susieR`, `jsonlite` and `openssl`. Supply the inputs documented in [data/raw/README.md](data/raw/README.md), or restore the private release according to its `PRIVATE_DATA.md`.

The tested Python 3.14 package versions are pinned in [requirements-analysis.lock.txt](requirements-analysis.lock.txt). Use that file instead of the minimum-version requirements when recreating the same software environment.

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

For regional SNP **and indel** pictures from an existing all-region run, add `--include-indels` to either regional plotter. The default retains the earlier SNP-only format. The separate [external fine-mapping workflow](docs/external_finemapping.md) extracts full regional summary statistics, acquires public reference LD, checks allele alignment, and fits exploratory SuSiE-RSS models. It does not repair unavailable cohort genotype calls or carrier associations.

## R code with editable options

Start with [R options and RStudio instructions](docs/R_options.md). Edit
[external fine-mapping settings](config/external_finemapping.json) or the
[R workflow example](config/R_workflow.example.json), then run:

```bash
Rscript code/R/finemap_external.R --config config/external_finemapping.json
```

The model supports joint SNP/indel, SNP-only or indel-only runs, configurable
credible-set coverage and model size, regional selection, workers, and
PNG/PDF/SVG exports. The [R source](code/R/finemap_external.R) and
[cohort R scripts](code/R/README.md) are included in both repositories.
Keep every real input and output in the private project.

## Validation

```bash
.venv/bin/python code/run_tests.py
Rscript code/R/test_analysis.R
Rscript code/R/test_audit.R
git diff --check
```

Tests use synthetic fixtures and do not establish results for a real cohort. The analysis preserves missing calls, validates phase, distinguishes genotype completeness from haplotype frequency, and reports unavailable estimates explicitly. Association and fine-mapping evidence do not by themselves establish causality.

See [AGENTS.md](AGENTS.md) for scientific and publishing requirements and [GitHub instructions](docs/github.md) for the separation between public code and private outputs.
