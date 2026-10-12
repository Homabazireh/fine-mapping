# Analysis results

Keep each analysis run in a new subdirectory. The configured cohort output is
`cohort_categories/`; the analysis requires its output directory to be new or
empty. No full cohort analysis has been run during report preparation.

`demo/` receives separate synthetic demonstration runs. `figures/` contains
report PNGs, PDFs and `figure_manifest.json` files in separate build directories
under `gwas/`, `demo/` or `cohort/`. The main analysis's optional `--plot` images
stay in `plots/` inside that run's output directory.

Build the local report from the project root:

```bash
.venv/bin/python code/build_report.py --render
```

The default report uses the configured real cohort run and your authoritative
category table. While sample ID aliases remain pending, it shows the GWAS and
category source only. `--demo` explicitly selects generated test participants;
older `figures/demo/build_*` folders are historical examples, not cohort results.
Use `--run-dir` to choose a different completed analysis.

`report.html` embeds the PNG images. PDF download links refer to separate files
under `figures/`. `report_content.qmd` is the generated include used by
`../report.qmd`; `report_context.json` records the selected source paths.

The renderer uses a temporary output directory and copies the HTML here after
successful rendering, preserving existing analysis results. Each report rebuild
gets a new figure directory. See `../docs/analysis.md` for output descriptions.
