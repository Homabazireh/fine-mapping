# Project settings

`settings.json` holds the active input/output paths and analysis defaults. JSON
keeps configuration readable with Python's standard library; no YAML package is
needed.

- `paths` defines `vcf`, `gwas`, `metadata`, `samples`, `out`, `demo_data` and
  `demo_results`. Relative values resolve from the project root, independent of
  the working directory.
- `analysis` contains the existing analysis parameters, including the
  `1e-8` P-value cutoff. The shipped defaults keep optional association tests
  and plotting disabled.

The main analysis accepts `--config path/to/settings.json`; a relative
configuration filename also resolves from the project root. A custom file can
provide only the settings being changed; missing settings use the built-in
defaults. Unknown keys or values with invalid types are rejected. Explicit
analysis CLI arguments override configuration values. Boolean settings can be
disabled explicitly with `--no-plot`, `--no-test-associations`,
`--no-allow-subset` or `--no-allow-filter-dot`.

The analysis still requires `--genome-build GRCh38` or `--genome-build GRCh37`
on the command line. This is a declaration that the inputs have already been
verified against the same build; configuration does not perform liftover.

`make_sample_template.py` uses the configured VCF but still requires an explicit
`--out`. Participant metadata is opt-in through `--metadata`; its presence in
settings does not automatically assign groups. The configured metadata source is your final
`data/raw/participant_groups_final.tsv`, with `name` identifiers and
`Diabetes_Group` categories. Prepare `data/intermediate/samples.tsv` using exact VCF IDs and independently verified aliases. Preserve source metadata and review records privately. If a saved preflight summary is used, verify that its metadata checksum matches the selected authoritative source.

## R workflow and external fine-mapping settings

- [external_finemapping.json](external_finemapping.json) controls the separate
  external R model and its figure export. Run
  `Rscript code/R/finemap_external.R --config config/external_finemapping.json`.
- [R_workflow.example.json](R_workflow.example.json) selects a native R analysis
  or plotting action. Inspect it using
  `Rscript code/R/run_configured_analysis.R --config config/R_workflow.example.json --dry-run`.

These two R configuration files resolve paths from the configuration directory;
`settings.json` retains its existing project-root path rule. CLI overrides use
the current working directory. Only the selected workflow action runs. Edit a
copy and use a new output directory for each model or display change.

[Full option table and RStudio examples](../docs/R_options.md) explain which
settings change model probabilities, which only change figures, and which
scientific quality requirements remain fixed.
