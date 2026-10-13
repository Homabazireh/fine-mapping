# Fine mapping

## Three-group comparison — 13 October 2026

The updated Python workflow and independent [R implementation](code/R/README.md)
compare Diabetes, Prediabetes and Control using shared marker panels, overall
three-group carrier tests, pairwise effects with 95% confidence intervals, and
separate Benjamini–Hochberg corrections for overall and pairwise tests.
See the [reproduction instructions](docs/haplotype_comparison.md),
[aggregate results](docs/figures/README.md), and
[agent instructions](AGENTS.md).

**The current data do not resolve which haplotypes differ.** All 20 selected
regions fail the requirements for usable haplotype panels. Consequently no
haplotype frequencies, group effects, p-values or causal variants are estimated.
The one published figure retains its correct meaning: genotype completeness.
Missing estimates are not evidence that the groups are equal. The evidence
audit explains why older exploratory fine-mapping results cannot resolve this.

The separate private repository `Homabazireh/fine-mapping-2` stores the project
data as private release archives. This public repository contains code and
reviewed aggregate results only.

This project finds regions highlighted by a GWAS, reads nearby variants from
your VCF, and compares local DNA patterns (**haplotypes**) between participant
groups. A haplotype is a set of variant alleles carried together on one
chromosome.

## Figures and private data

All saved figures and complete analysis results are available in the private
[Fine mapping 2 repository](https://github.com/Homabazireh/fine-mapping-2).
Its `all_figures/` catalog provides browsable copies, and its private release
contains all project inputs and detailed results with restoration instructions.

This public repository retains the full Python and R code, documentation and
[reviewed aggregate result tables](docs/figures/README.md). Figure files and
reports containing embedded figures have been removed from the current public
branch at the user's request. Do not republish them here.

## Start from a fresh checkout

Use Python 3.10 or newer. Run these commands from the repository root:

```bash
python3 -m venv .venv
.venv/bin/python code/run_tests.py
.venv/bin/python code/run_demo.py
```

The core tests and synthetic demo need no private inputs or third-party Python
packages. Tests requiring local GWAS data or optional plotting packages are
skipped when those inputs are unavailable. The demo creates generated inputs
under `data/intermediate/demo/` and a completed run under `results/demo/`.

The repository contains code, settings, report source, literature notes and the
reviewed aggregate tables in `docs/figures/`. Raw data, participant
tables, review files, analysis results, rendered reports, archived originals and
local tool installations are ignored by Git and are **not included in a clone**.

## Current local cohort status

In the prepared analysis workspace, `results/report.html` contains figures and
explanations for:

- An overview of your real `EUR_diagram_t2d_1e-8.tsv` GWAS file.
- Participant counts displayed as **Diabetes, Prediabetes and Control**.
- Lead-genotype completeness for **all 147 participants**, using 136 exact
  metadata–VCF ID matches and 11 user-confirmed aliases.
- The completed cohort analysis and why none of its 20 selected regions
  produced a haplotype estimate.

Your final supplied table, `data/raw/participant_groups_final.tsv`, is the
authoritative category source: **49 Diabetes, 48 Prediabetes and 50 Control**.
Its columns are `External Name`, `Diabetes_Group` and `name`; the code reads
groups from `Diabetes_Group` and participant IDs from `name`. Source values are
`T2D`, `Prediabetes` and `Control`, with `T2D` displayed as Diabetes. Compared
with the previous table, 92 participants change phenotype group and another
20 have a label-only change from `NORM` to `Control`. All 11 proposed identifier
pairs are now confirmed and recorded in
`data/intermediate/participant_review/verified_id_map.tsv`.

The group-QC figure uses all 49 Diabetes, 48 Prediabetes and 50 Control
participants. These comprise 45/44/47 exact matches plus 4/4/3 verified aliases,
respectively. No participant identifiers remain unresolved.

An indexed check of the 20 default lead variants also found that none meets
the current 90% genotype-completeness requirement. The highest is 126/147
(85.71%). Missing genotypes remain missing; the report records this limitation.

The full analysis completed in `results/cohort_categories/` with all 147
participants and no exclusions. It produced **no haplotype frequency estimates
or association tests**: 19 regions were skipped for insufficient lead-genotype
completeness and one for unconfirmed palindromic lead-allele orientation. The
completed run resolves the grouping question; the available genotypes still
prevent the intended haplotype comparison.

## Why an older figure says “synthetic demo”

`results/demo/` contains runs on artificial participants generated by
`code/run_demo.py`. Haplotype and phase figures under
`results/figures/demo/build_*/` illustrate those generated genotypes. Their
synthetic labels are intentional. A demo build can also contain a real GWAS or
participant-count overview; each figure's `scope` in its manifest identifies
its source.

An old editor tab can still point to a previous demo build. Open the current
`results/report.html` and check `results/report_context.json`: its
`figure_manifest` field points to the active figure set. The current report has
`synthetic: false`; its `run_dir` selects `results/cohort_categories/`.
Genotype completeness and the completed cohort outcome appear in separate
sections.

Full reports, working figure directories and participant data stay local and
are ignored by Git. All figure copies are shared only through the private `fine-mapping-2` repository.

## Build figures and the report

Supply the local GWAS and participant table at the paths in
`config/settings.json`, or change those paths to your own files. The report
builder always reads the configured GWAS, including when `--demo` is selected.
Participant-category figures require the configured metadata table.

Install Matplotlib and [Quarto](https://quarto.org/docs/download/index.html),
then run:

```bash
.venv/bin/python -m pip install matplotlib
.venv/bin/python code/build_report.py \
  --preflight data/intermediate/confirmed_groups_qc/summary.json --render
```

This uses the configured cohort run when complete, plus the GWAS,
authoritative participant table and selected genotype-completeness check. It
saves `results/report.html` and a new folder under `results/figures/`.

Use `--preflight path/to/summary.json` to choose a fresh genotype check that
uses the same participant table. The current check uses exact VCF ID matches
and the 11 aliases supplied through `--id-map`. The report records both match
types and uses their combined counts for genotype quality control.

To include a synthetic example alongside the configured local GWAS:

```bash
.venv/bin/python code/run_demo.py
.venv/bin/python code/build_report.py --demo \
  --preflight data/intermediate/confirmed_groups_qc/summary.json --render
```

## Which files should I use?

| File or folder | What it is for |
| --- | --- |
| `code/compare_haplotypes.py` | Main analysis: select regions, read genotypes, compare haplotypes |
| `code/make_sample_template.py` | Prepare a participant table for checking |
| `code/check_lead_genotypes.py` | Check genotype completeness at the selected GWAS leads |
| `code/make_figures.py` | Draw the GWAS and haplotype figures |
| `code/comparison_figures.py` | Draw paper-inspired cohort and regional panels from saved analysis tables |
| `code/export_panel_c.py` | Export only panel c at 600 dpi, with vector copies and a focused report |
| `code/build_report.py` | Put explanations and figures into the report |
| `code/run_demo.py` | Run the analysis on synthetic data |
| `config/settings.json` | Change input paths and analysis settings |
| `data/raw/` | Original GWAS, participant metadata and links to the VCF |
| `data/intermediate/` | Participant checks and synthetic inputs |
| `results/` | Analysis tables, figures and rendered report |
| `report.qmd` | Quarto report title, layout and included content |
| `docs/analysis.md` | Plain-language workflow and detailed commands |
| [literature_review/](literature_review/) | Research papers, references and reading notes |
| `archive/` | Local earlier copies, including the old `(1)` files; not distributed |

Open the scripts under **`code/`** for editing. In the original workspace,
`data/raw/full_LR-WGS.vcf.gz` links to the VCF in the parent directory.
On another machine, supply that input or set `paths.vcf` to its location.

## Use the real participants

In the prepared workspace, `data/intermediate/samples.tsv` contains all 147
verified assignments with columns `sample` and `group`. The 11 confirmed
aliases are recorded in `data/intermediate/participant_review/verified_id_map.tsv`;
all sample-table identifiers match the VCF. These local files are not distributed
with the repository. The [analysis guide](docs/analysis.md#participant-groups)
explains how to prepare a new table from the supplied metadata and verified map.

After checking the groups and confirming both inputs use GRCh38:

```bash
.venv/bin/python code/compare_haplotypes.py \
  --samples data/intermediate/samples.tsv \
  --out results/cohort_run_01 \
  --genome-build GRCh38

.venv/bin/python code/build_report.py \
  --run-dir results/cohort_run_01 \
  --preflight data/intermediate/confirmed_groups_qc/summary.json --render
```

Choose a new output folder for each analysis. The configured GWAS cutoff is
`P ≤ 1e-8`. Change defaults in `config/settings.json` or use command-line options.

## Check the code

```bash
.venv/bin/python code/run_tests.py
```

Figures need Matplotlib; rendering needs Quarto. The original workspace has
local installations, which are not distributed. See the
[analysis guide](docs/analysis.md) for setup and optional statistical tests.
