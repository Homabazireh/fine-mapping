# Analysis results

Keep each analysis run in a new subdirectory. The configured cohort output is
`cohort_categories/`; the analysis requires its output directory to be new or
empty. No full cohort analysis has been run during report preparation.

Open `report.html` for the current local report. `report_context.json` records
the selected inputs; follow its `figure_manifest` field to the active PNGs,
PDFs and captions. Each rebuild creates a new figure directory, so a previously
opened editor tab can still show an older build.

The current report has `synthetic: false` and no completed cohort `run_dir`.
It contains three real-data figure types:

- The supplied public GWAS overview.
- Source-table counts: 49 Diabetes, 47 Prediabetes and 51 Control.
- Lead-genotype completeness for 46 Diabetes, 43 Prediabetes and 47 Control
  participants whose metadata IDs exactly match the VCF.

The 11 other VCF identifiers remain unassigned. No real cohort haplotype
frequency or phase results are available yet. Genotype-completeness percentages
are quality-control results, not haplotype frequencies.

`demo/` contains runs on artificial participants and genotypes made by
`code/run_demo.py`. Haplotype and phase plots in `figures/demo/build_*/` describe
those synthetic examples. A demo report can also include real GWAS and
source-category summaries; the manifest's per-figure `scope` distinguishes
these inputs. A “SYNTHETIC DEMO” label on an old example plot is intentional.

Other figure builds live under `figures/gwas/` or `figures/cohort/`. The main
analysis's optional `--plot` images stay in `plots/` inside its run directory.

Build the local report from the project root:

```bash
.venv/bin/python code/build_report.py --render
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
