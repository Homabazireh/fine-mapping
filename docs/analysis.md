# Understand and run the analysis

The question is: **Do participant groups carry different local DNA patterns
near variants associated with type 2 diabetes?**

The public GWAS chooses the regions. The VCF supplies participant genotypes.
The analysis counts haplotypes in each group, and the report explains the
figures. The prepared local report uses the real GWAS and the authoritative Category table.
The real cohort comparison awaits confirmation of 11 VCF identifier aliases.

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
| Supplied participant table | `data/raw/participant_categories.tsv` |
| Final, verified sample groups | `data/intermediate/samples.tsv` — pending |

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

Your authoritative table has **147 people: 51 NORM, 47 Prediabetes and 49
T2D**. Its categories supersede the 88 differences in the older ID lists.
There are **136 exact VCF matches** and **11 proposed aliases awaiting
confirmation**. Review the remaining ID differences before the real run:

- `sample_comparison.tsv`: compare both group sources for each VCF sample.
- `metadata_ids_not_in_vcf.tsv`: check the 11 unmatched identifiers.
- `summary.json`: see aggregate counts.

These files are under `data/intermediate/participant_review/` in the prepared
workspace and are not included in a clone.

Prepare a table from the supplied metadata:

```bash
.venv/bin/python code/make_sample_template.py \
  --metadata data/raw/participant_categories.tsv \
  --out data/intermediate/samples_for_review.tsv
```

Exact matches receive the authoritative categories. Unmatched VCF IDs receive
`UNRESOLVED`. Verify those ID mappings, then save the completed assignments as
`data/intermediate/samples.tsv`. It needs `sample` and `group`
columns, exact VCF IDs, and at least two groups. Category names are preserved;
the helper refuses to overwrite existing output files.

For independently verified aliases, provide a one-to-one TSV with
`metadata_sample` and `vcf_sample` columns using
`--id-map data/intermediate/verified_id_map.tsv`. The code never guesses aliases.
Omit `--metadata` to start with all groups unassigned.

## Current genotype completeness

An indexed check of the 20 default lead variants found an allele-matching PASS
record for each, but none had usable genotypes for at least 90% of the 147 VCF
participants at DP ≥ 10 and GQ ≥ 20. The highest was 126/147 (85.71%). No
reference-homozygous genotypes were present at those lead sites. Missing calls
are retained as missing.

Reproduce this aggregate check, which does not need phenotype assignments:

```bash
.venv/bin/python code/check_lead_genotypes.py \
  --out data/intermediate/cohort_preflight_repeat
```

See `data/intermediate/cohort_preflight/lead_genotype_qc.tsv` and `summary.json`.
The checker needs `bcftools`, the local VCF and its index. Install or provide
them on a new machine. Choose a new output folder for each repeat so previous
checks are preserved.
Resolving the participant IDs alone will therefore not produce haplotype
estimates under the current settings. The VCF's genotype completeness needs
investigation before interpreting group comparisons.

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

## Figures and Quarto report

Report generation requires Matplotlib and the configured local GWAS file.
Participant-category figures also require the configured metadata table;
cohort figures require a completed local analysis. HTML rendering additionally
requires Quarto. None of the local data, output files or tool installations is
downloaded with the repository.

Rebuild the current report from the configured real cohort run (or the GWAS
and authoritative category source while cohort analysis is pending):

```bash
.venv/bin/python code/build_report.py --render
```

Make a report from a completed real run:

```bash
.venv/bin/python code/build_report.py \
  --run-dir results/cohort_run_01 --render
```

Open `results/report.html`. Read the figures as follows:

| Figure | How to read it |
| --- | --- |
| GWAS overview | Each dot is a supplied variant. Higher points have smaller P-values. The lower panel zooms in on the threshold. |
| Participant categories | Counts from your authoritative table: T2D, Prediabetes and NORM. These are source counts, before VCF ID matching. |
| Participants with resolved haplotypes | Each cell shows people with usable haplotype pairs, out of their whole group. |
| Haplotype frequencies | Each cell shows the fraction of resolved chromosome copies carrying that pattern. Zero means absent among resolved copies; `NA` means no estimate. |

Read frequencies alongside phase resolution: different rates of unresolved
participants can affect group comparisons. Synthetic figures illustrate the
software and do not describe your cohort.

`report.qmd` controls layout. `code/build_report.py` creates explanations and
figure references in `results/report_content.qmd`. Edit the builder for wording
that should persist after a rebuild. `code/make_figures.py` draws the figures.

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
