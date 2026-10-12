# Understand and run the analysis

The question is: **Do participant groups carry different local DNA patterns
near variants associated with type 2 diabetes?**

The public GWAS chooses the regions. The VCF supplies participant genotypes.
The analysis counts haplotypes in each group, and the report explains the
figures. The prepared local report uses the real GWAS and the final supplied
participant table, with groups read from its `Diabetes_Group` column.
All 147 participant identities are matched: 136 by exact identifiers and 11
through user-confirmed aliases. Genotype completeness determines which regions
can produce haplotype estimates. The cohort analysis has completed; no selected
region produced an eligible haplotype panel.

Run the commands below from the repository root (`Fine_mapping/new_analysis/`
in the original workspace, or your `fine-mapping/` checkout).

Raw inputs, participant tables, generated outputs and local environments are
ignored by Git. A fresh checkout includes the code and report source, but not
`results/report.html` or the local cohort review files. The synthetic demo and
core tests work without those data; see [setup, demo and tests](#setup-demo-and-tests).

## Five steps in the code

Read `code/compare_haplotypes.py` starting at `run()`. It coordinates these steps:

1. **Check the inputs.** Read the VCF header and participant groups. IDs must
   match exactly, and every included participant needs a known group.
2. **Choose regions from the GWAS.** Keep variants with `P ≤ 1e-8`. Group nearby
   variants and choose the lowest-P variant in each group as its lead.
3. **Read nearby genotypes.** Scan the VCF once and retain eligible SNPs and
   small insertions/deletions around each lead.
4. **Build and count haplotypes.** Select markers, identify each participant's
   two haplotypes where possible, and count their copies in each group.
5. **Save the results.** Write the selected markers, frequencies, unresolved
   calls, and run settings. Separate scripts draw figures and build the report.

The workflow compares local haplotypes. Statistical fine-mapping, liftover,
and filling in missing genotypes are outside its scope.

## Useful terms

| Term | Meaning here |
| --- | --- |
| GWAS | A study testing associations between variants and a trait |
| VCF | A file containing variants and participant genotypes |
| Lead variant | The variant with the smallest P-value in a nearby group |
| Locus | A selected region of the genome |
| Marker | A variant used to describe a local DNA pattern |
| Haplotype | Marker alleles carried together on one chromosome |
| Phase | Information about which alleles occur on the same chromosome |
| Resolved participant | Someone whose two haplotypes can be determined at the selected markers |
| Haplotype frequency | Copies of a haplotype divided by all resolved chromosome copies in that group |

## Inputs and settings

| Input | Location |
| --- | --- |
| Filtered GWAS used by default | `data/raw/EUR_diagram_t2d_1e-8.tsv` |
| Alternative full GWAS | `data/raw/EUR_diagram_t2d.hg38.tsv` |
| Cohort VCF | `data/raw/full_LR-WGS.vcf.gz` |
| Final supplied participant table | `data/raw/participant_groups_final.tsv` |
| Final, verified sample groups | `data/intermediate/samples.tsv` |
| Confirmed identifier aliases | `data/intermediate/participant_review/verified_id_map.tsv` |

In the original workspace, the VCF path is a symbolic link to
`full_LR-WGS.vcf.gz` in the parent directory. That link and its target are not
distributed. Supply your own VCF and set `paths.vcf` in `config/settings.json`.
The main analysis reads it sequentially and does not use the `.csi` index.

The filtered GWAS contains **64,929 rows**: **54,060 with P ≤ 1e-8** and
**10,869 above this cutoff**, up to `5e-8`. The code applies the cutoff again.
It preserves very small P-values such as `3.35e-1003` using `Decimal`.

The main defaults in `config/settings.json` are:

| Setting | Default | Meaning |
| --- | --- | --- |
| `p_threshold` | `1e-8` | Keep P-values at or below this value |
| `lead_gap` | `150000` | Maximum gap between consecutive variants in a cluster, in bases |
| `top_loci` | `20` | Analyze the 20 clusters with the strongest lead signals |
| `flank` | `50000` | Read 50,000 bases on each side of a lead |
| `max_markers` | `8` | Use up to eight markers per locus, including the lead |
| `min_dp` / `min_gq` | `10` / `20` | Minimum genotype read depth and quality |
| `min_call_rate` | `0.9` | Require usable genotypes in at least 90% of included participants |

Command-line options override settings. Paths in settings are project-relative;
explicit relative command-line paths use your working directory.
See [configuration notes](../config/README.md) for complete settings rules.

## Participant groups

Your final authoritative table has **147 people: 49 T2D, 48 Prediabetes and 50
Control**. The columns are `External Name`, `Diabetes_Group` and `name`. The
code reads participant IDs from `name` and assignments from `Diabetes_Group`;
`External Name` is retained as supplied metadata. This final table supersedes
the earlier participant table and ID lists. Compared with the previous table,
**92 participants change phenotype group** and **20 additional labels change
from NORM to Control** without changing the displayed phenotype group.

Figures and report tables display `T2D` as **Diabetes**, `Prediabetes` as
**Prediabetes**, and `Control` as **Control**. Input category values are preserved.
The report shows both the source category and display label, with exact VCF
matches and verified aliases counted separately for each group.

There are **136 exact VCF matches**: 45 Diabetes, 44 Prediabetes and 47 Control.
The remaining **11 aliases are now confirmed**: 4 Diabetes, 4 Prediabetes and
3 Control. Thus all 49 Diabetes, 48 Prediabetes and 50 Control participants
are assigned. The review records include:

- `sample_comparison.tsv`: compare group sources for each VCF sample.
- `category_assignment_change_audit.tsv`: record changes from the previous table.
- `metadata_ids_not_in_vcf.tsv`: retain identifiers that lacked an exact match.
- `verified_id_map.tsv`: record the 11 confirmed metadata-to-VCF aliases.
- `summary.json`: see aggregate counts.

These files are under `data/intermediate/participant_review/` in the prepared
workspace and are not included in a clone.

The configured `data/intermediate/samples.tsv` contains all 147 assignments.
To prepare a new table from the supplied metadata and verified map:

```bash
.venv/bin/python code/make_sample_template.py \
  --metadata data/raw/participant_groups_final.tsv \
  --id-map data/intermediate/participant_review/verified_id_map.tsv \
  --out data/intermediate/samples_for_review.tsv
```

Exact matches and verified aliases receive the authoritative categories; the
output identifiers are the VCF IDs. Any unmatched VCF IDs in a future table
receive `UNRESOLVED`. Final analysis input needs `sample` and `group` columns,
verified VCF IDs, and at least two groups. Category names are preserved; the
helper refuses to overwrite existing output files.

For independently verified aliases, provide a one-to-one TSV with
`metadata_sample` and `vcf_sample` columns using
`--id-map data/intermediate/participant_review/verified_id_map.tsv`. The code
never guesses aliases.
Omit `--metadata` to start with all groups unassigned.

## Current genotype completeness

An indexed check of the 20 default lead variants found an allele-matching PASS
record for each, but none had usable genotypes for at least 90% of the 147 VCF
participants at DP ≥ 10 and GQ ≥ 20. The highest was 126/147 (85.71%). No
reference-homozygous genotypes were present at those lead sites. Missing calls
are retained as missing.

Create the check and group summaries using the final supplied table:

```bash
.venv/bin/python code/check_lead_genotypes.py \
  --metadata data/raw/participant_groups_final.tsv \
  --id-map data/intermediate/participant_review/verified_id_map.tsv \
  --out data/intermediate/confirmed_groups_qc
```

See `data/intermediate/confirmed_groups_qc/lead_genotype_qc.tsv` and `summary.json`.
The checker needs `bcftools`, the local VCF and its index. Install or provide
them on a new machine. Choose a new output folder for each repeat so previous
checks are preserved.
The now-confirmed participant IDs alone do not ensure haplotype estimates
under the current settings. The VCF's genotype completeness needs
investigation before interpreting group comparisons.

With `--metadata` and the verified `--id-map`, the checker also reports
lead-genotype completeness for all 49 Diabetes, 48 Prediabetes and 50 Control
participants. It records exact counts of 45/44/47 and verified-alias counts of
4/4/3 separately, with their combined counts as the QC denominators. All 147
participants are assigned. The check measures genotype completeness and does
not calculate haplotype frequencies.

## Completed cohort result

`results/cohort_categories/` contains the completed analysis with all 147
participants included and none excluded. It scanned 33,244,057 VCF records,
including 16,110 records within the selected windows. Of the 20 selected
regions, 19 were skipped because the lead genotype failed the completeness
threshold; one was skipped because its palindromic lead-allele orientation
was unconfirmed. The separate QC check shows that all 20 leads fall below
90% usable-genotype completeness.

No haplotype panels, frequency estimates or carrier association tests were
produced. No high-quality reference-homozygous genotypes were observed in the
selected windows. Missing calls were retained as missing. The sample mapping
is complete; these results identify genotype and allele-orientation limits
on the requested comparison, and do not establish equal group frequencies.
See `lead_loci.tsv`, `run_summary.json` and `COMPLETED.txt` in the run directory.

## Run the real comparison

After checking the groups and confirming both inputs use GRCh38:

```bash
.venv/bin/python code/compare_haplotypes.py \
  --samples data/intermediate/samples.tsv \
  --out results/cohort_run_01 \
  --genome-build GRCh38
```

Each run requires a new or empty output folder. Use `--top-loci 0` for all
clusters, or `--help` to see all options. Groups may have different sizes.
Every VCF sample needs an assignment unless you deliberately use `--allow-subset`.

The analysis runs directly from the VCF, GWAS and sample table. It does not
call `check_lead_genotypes.py` or read its saved summaries. Omitting that
separate check leaves the analysis's DP, GQ, call-rate and phase requirements
active. Changing those thresholds is a separate analysis decision.

## Figures and Quarto report

### Standalone panel c

The [published panel c](figures/README.md) shows lead-genotype completeness
across the 20 selected regions. It has a 600-dpi PNG, vector PDF/SVG and a
self-contained one-figure HTML report. Reproduce it from a completed local run:

```bash
.venv/bin/python code/export_panel_c.py \
  --run-dir results/cohort_no_preflight \
  --out results/cohort_no_preflight/panel_c_new
```

Use a new or empty output folder. This export reads saved aggregate tables;
it does not run or load a separate VCF precheck. The original analysis's
filters remain unchanged. The current output is at
`results/cohort_no_preflight/panel_c/`.

### Full report

Report generation requires Matplotlib and the configured local GWAS file.
Participant-category figures also require the configured metadata table;
cohort figures require a completed local analysis. HTML rendering additionally
requires Quarto. None of the local data, output files or tool installations is
downloaded with the repository.

Rebuild the current report from the configured real cohort run, GWAS,
authoritative category source and confirmed-participant QC:

```bash
.venv/bin/python code/build_report.py \
  --preflight data/intermediate/confirmed_groups_qc/summary.json --render
```

To save a report beside a completed direct run, without loading any saved
preflight summary, use:

```bash
.venv/bin/python code/build_report.py \
  --run-dir results/cohort_no_preflight \
  --no-preflight --output-dir results/cohort_no_preflight --render
```

This puts `report.html`, its source and figure exports in the specified folder.
It preserves the default report and the run's analysis tables. The report
records that the standalone preflight was omitted; internal analysis filters
remain as recorded in `run_summary.json`.

If you repeat the check under a new output directory, supply its summary:

```bash
.venv/bin/python code/build_report.py \
  --preflight data/intermediate/cohort_preflight_repeat/summary.json --render
```

Choose the summary explicitly so it matches the final participant source.
Without `--preflight`, the builder can use a saved default summary from an
older check. A supplied missing summary is an error, and group-QC metadata
checksums must match the configured participant table.

Make a report from a completed real run:

```bash
.venv/bin/python code/build_report.py \
  --run-dir results/cohort_run_01 \
  --preflight data/intermediate/confirmed_groups_qc/summary.json --render
```

Open `results/report.html`. Read the figures as follows:

| Figure | How to read it |
| --- | --- |
| GWAS overview | Each dot is a supplied variant. Higher points have smaller P-values. The lower panel zooms in on the threshold. |
| Participant categories | Diabetes, Prediabetes and Control counts from the final table; the source values are T2D, Prediabetes and Control. These are totals before VCF ID matching. |
| Lead-genotype completeness by group | Usable lead genotypes divided by all matched participants in each group: 49 Diabetes, 48 Prediabetes and 50 Control. Exact matches and verified aliases contribute to these denominators and are reported separately. These are QC percentages, not haplotype frequencies. |
| Cohort comparison overview | Lettered panels show the analysis workflow, included participants and usable genotypes at all selected leads. Black diamonds show pooled completeness; the dashed line marks the pooled requirement. Values come from saved main-analysis tables, so no separate precheck is required. |
| Regional marker comparison | Group rows share a physical position axis. Left: called participants / group size. Right: ALT copies / called chromosomes. Sites with no calls have zero completeness and no allele-frequency estimate. The example is the unique exact biallelic SNP lead with highest pooled coverage, with ties resolved by GWAS lead-table order. These points do not represent phase, shared ancestry or tested differences. |
| Participants with resolved haplotypes | Each cell shows people with usable haplotype pairs, out of their whole group. |
| Haplotype frequencies | Each cell shows the fraction of resolved chromosome copies carrying that pattern. Zero means absent among resolved copies; `NA` means no estimate. |

Read frequencies alongside phase resolution: different rates of unresolved
participants can affect group comparisons. Synthetic figures illustrate the
software and do not describe your cohort.

`report.qmd` controls layout. `code/build_report.py` creates explanations and
figure references in `results/report_content.qmd`. Edit the builder for wording
that should persist after a rebuild. `code/make_figures.py` draws the figures;
`code/comparison_figures.py` prepares the paper-inspired cohort and regional
panels. These use the layout concepts in [Figueroa et al. (2024), Figure 1](https://www.nature.com/articles/s41588-024-01719-5/figures/1)
with this project's own data. The accompanying `*_plot_data.tsv` files expose
the plotted numerators and denominators. No supplied or generated chromosome
patterns are substituted for unresolved haplotypes.

Every build gets a new folder under `results/figures/` with PNGs, PDFs and a
`figure_manifest.json` recording sources and captions. The latest selection is
in `results/report_context.json`. The HTML embeds its PNG images; keep its
corresponding `figures/` folder alongside it for PDF links to work.

`--demo` selects the latest completed synthetic run. Without either run option,
the report uses the configured completed cohort run, or shows the GWAS and
participant review if none is available. An explicitly requested missing run
is an error. Rendering uses a temporary directory to preserve existing results;
it never starts a cohort analysis.

Even with `--demo`, the public-GWAS overview reads the configured local GWAS;
the synthetic run supplies only the example haplotype results. Running the
synthetic analysis itself does not require that local GWAS.

## Main output tables

| File | What to look for |
| --- | --- |
| `lead_loci.tsv` | Selected regions and reasons for any skipped regions |
| `selected_markers.tsv` | Variants used to define each haplotype |
| `sample_haplotypes.tsv` | Each participant's haplotype pair or unresolved status |
| `haplotype_group_frequencies.tsv` | Copy counts, frequencies and denominators |
| `group_phase_qc.tsv` | Resolved participants and exclusion reasons by group |
| `haplotype_comparisons.tsv` | Group comparisons and any requested tests |
| `run_summary.json` | Input files, settings, counts and warnings |

Runs also save `samples_used.tsv`, `selected_marker_genotypes.tsv.gz`, variant
quality tables, and `regions.bed`. TSV positions are 1-based; BED intervals are
0-based and exclude the end position. Keep a run's files together.

## Technical details when needed

- GWAS columns include `Chromsome` (the supplied spelling), `Position`,
  `EffectAllele`, `NonEffectAllele` and `Pval`. `Beta` and `SE` break P-value
  ties; `liftover_strand` helps check allele orientation.
- Nearby markers come from the VCF and may be absent from the filtered GWAS.
  Select the full GWAS with `--gwas data/raw/EUR_diagram_t2d.hg38.tsv` if needed.
  `--p-threshold 5e-8` admits supplied rows up to that P-value; other checks apply.
- Site `FILTER=PASS` is required by default. A/T and C/G lead SNPs and
  negative-strand liftover leads need verified forward-reference alleles before
  using `--gwas-alleles-forward-confirmed`.
- `--genome-build` records your declaration. Both inputs must already have
  matching builds and compatible allele representations.
- `--test-associations` requests exploratory carrier tests. It requires SciPy
  and verified reference-genotype callability, declared through
  `--reference-genotypes-confirmed`. These tests do not adjust for covariates.
- `--plot` saves optional marker-pattern PNGs inside a run. Report frequency and
  phase figures are generated separately by `make_figures.py`.

## Setup, demo and tests

From a fresh checkout, use Python 3.10 or newer. The core analysis, synthetic
demo and core tests use only the Python standard library:

```bash
python3 -m venv .venv
.venv/bin/python code/run_tests.py
.venv/bin/python code/run_demo.py
```

Tests depending on unavailable local GWAS data or optional plotting packages
are skipped. To enable figures and plotting tests:

```bash
.venv/bin/python -m pip install matplotlib
```

Install [Quarto](https://quarto.org/docs/download/index.html) to render HTML.
The builder finds it on `PATH` or at `.tools/quarto/bin/quarto`.
`requirements.txt` also lists SciPy for optional association tests. The original
workspace's `.venv` and `.tools/` are local installations and are not part of
the repository.

The demo writes synthetic inputs to `data/intermediate/demo/` and a fresh run
to `results/demo/`. Its two positional arguments can override those folders.
Use `--config path/to/settings.json` to load another settings file.
