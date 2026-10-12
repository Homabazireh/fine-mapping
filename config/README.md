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
settings does not automatically assign groups. The configured metadata source is your authoritative
`data/raw/participant_categories.tsv`. The final `data/intermediate/samples.tsv`
awaits confirmation of 11 VCF ID aliases; the category-source choice is settled.
