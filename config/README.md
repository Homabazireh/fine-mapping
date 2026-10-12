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
`Diabetes_Group` categories. The final `data/intermediate/samples.tsv` contains
147 assignments, using 136 exact ID matches and 11 user-confirmed aliases from
`data/intermediate/participant_review/verified_id_map.tsv`.

Rebuild the current report with
`--preflight data/intermediate/confirmed_groups_qc/summary.json` to select the
genotype check computed from that final table. Earlier summaries may refer to
superseded metadata and must not supply the current group denominators.
