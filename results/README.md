# Analysis results

Keep each analysis run in a new subdirectory. The configured cohort output is
`cohort_categories/`; the analysis requires its output directory to be new or
empty. The full cohort analysis has completed in this directory with all 147
participants included and none excluded.

Open `report.html` for the current local report. `report_context.json` records
the selected inputs; follow its `figure_manifest` field to the active PNGs,
PDFs and captions. Each rebuild creates a new figure directory, so a previously
opened editor tab can still show an older build.

The current report has `synthetic: false` and selects `cohort_categories/` as
its completed cohort `run_dir`. It includes five real-data figure types:

- The supplied public GWAS overview.
- Source-table counts: 49 Diabetes, 48 Prediabetes and 50 Control from the final
  `data/raw/participant_groups_final.tsv` table's `Diabetes_Group` column.
- Lead-genotype completeness for all 49 Diabetes, 48 Prediabetes and 50 Control
  participants, using 136 exact matches and 11 user-confirmed aliases.
- A paper-inspired overview combining the workflow, included participants and
  all 20 lead sites, calculated from the completed analysis tables.
- A regional comparison at L0010, selected by highest pooled lead completeness:
  528 biallelic SNPs on a shared physical position axis, with genotype completeness
  and observed ALT-copy fractions in separate panels. These are recorded-call
  summaries; they do not represent phased haplotypes or unbiased population frequencies.

Both added figures use the panel organization and regional-comparison idea in
[Figueroa et al. (2024), Figure 1](https://www.nature.com/articles/s41588-024-01719-5/figures/1).
They are also included in `cohort_no_preflight/report.html`, without loading a
separate VCF precheck. A convenient export bundle is saved at
`cohort_no_preflight/paper_inspired_figures/`, with PNGs, PDFs and the exact plot
tables. The generated figure manifest records the example-selection rule and
source files. No new genotype or haplotype analysis was run to create these plots.

No participant identifiers remain unresolved. Exact matches number 45/44/47
and verified aliases 4/4/3 for Diabetes/Prediabetes/Control, respectively.
Genotype-completeness percentages are quality-control results, not haplotype
frequencies. Read the cohort section for the selected analysis status.
The final table changes 92 phenotype assignments relative to the previous
table, plus 20 label-only changes from `NORM` to `Control`.

The completed run produced no eligible haplotype panels, frequency estimates
or carrier association tests. Of 20 selected regions, 19 were skipped for
insufficient lead-genotype completeness and one for unconfirmed palindromic
lead-allele orientation. The separate QC check shows that all 20 leads fail
the 90% completeness threshold. The run scanned 33,244,057 VCF records,
including 16,110 within the selected windows. See `cohort_categories/lead_loci.tsv`
and `run_summary.json` for the reasons and settings; `COMPLETED.txt` records
successful completion. These limitations concern the genotype data and allele
orientation; all participant identifiers and categories are assigned.

`demo/` contains runs on artificial participants and genotypes made by
`code/run_demo.py`. Haplotype and phase plots in `figures/demo/build_*/` describe
those synthetic examples. A demo report can also include real GWAS and
source-category summaries; the manifest's per-figure `scope` distinguishes
these inputs. A “SYNTHETIC DEMO” label on an old example plot is intentional.

Other figure builds live under `figures/gwas/` or `figures/cohort/`. The main
analysis's optional `--plot` images stay in `plots/` inside its run directory.

Build the local report from the project root:

```bash
.venv/bin/python code/build_report.py \
  --preflight data/intermediate/confirmed_groups_qc/summary.json --render
```

The default report uses the configured real cohort run and your authoritative
category table. When no completed cohort run is available, it shows the GWAS,
category source and any selected genotype check. Use `--preflight path/to/summary.json`
to select a refreshed check. `--demo` explicitly selects generated test participants.
Use `--run-dir` to choose a different completed analysis.

`report.html` embeds the PNG images. PDF download links refer to separate files
under `figures/`. `report_content.qmd` is the generated include used by
`../report.qmd`; `report_context.json` records the selected source paths.

The renderer uses a temporary output directory and copies the HTML here after
successful rendering, preserving existing analysis results. Each report rebuild
gets a new figure directory. See `../docs/analysis.md` for output descriptions.

Generated results, figures, rendered reports and participant data are local
artifacts ignored by Git. They are not included when the source code is
published to GitHub. This README is the documentation exception.

Selected aggregate PNG/PDF copies are also published in
[`docs/figures/`](../docs/figures/README.md), with captions and a provenance
manifest. This does not publish the underlying participant or genotype tables.
