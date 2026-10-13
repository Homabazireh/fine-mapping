# Working on this analysis

## Scientific requirements

- Use `config/settings.json` and the final participant table as the input sources. The current verified cohort contains 49 Diabetes, 48 Prediabetes and 50 Control participants; do not reuse superseded category lists.
- Preserve missing genotypes as missing. Never replace absent calls with reference alleles. Do not relax genotype, phase, allele-orientation or reference-provenance requirements merely to obtain results.
- Use one common, phenotype-blind marker panel per locus. Treat participants, rather than their two correlated chromosome copies, as the independent units in carrier association tests.
- Report overall and pairwise test families separately, with multiple-testing correction, effect direction and uncertainty. Distinguish descriptive chromosome-copy frequencies from carrier association tests.
- Report unavailable estimates as unavailable, with reasons. No eligible haplotypes does not mean equal frequencies, no association, or no causal variants.
- Distinguish external GWAS evidence, exploratory fine-mapping and cohort associations. Never promote unvalidated legacy PIPs to causal evidence.
- This requested publication contains only panel C. PNG, PDF and SVG are formats of the same figure; do not add other figures to its report.

## Reproducibility

- Keep the Python and independent R implementations consistent. Record package versions and input provenance. Preserve previous runs; use a new output directory.
- Run `.venv/bin/python code/run_tests.py`, the R tests documented in `code/R/README.md`, and `git diff --check` after relevant changes.
- Add meaningful synthetic regression tests for statistical, genotype and phase behavior. Never present synthetic results as cohort results.

## Publishing and data

- `Homabazireh/fine-mapping` is public: publish code, documentation and reviewed aggregate results only. Raw genotypes, participant tables, aliases and credentials must remain excluded.
- The user authorized a separate private repository, `Homabazireh/fine-mapping-2`, for project data. Verify that destination remains private before any data upload. Do not add collaborators or change visibility.
- Large project data belong in checked, checksum-verified private release archives. Keep credentials, Git internals, installed environments and caches out of archives. Include a file manifest and restoration instructions.
- Inspect staged files and the destination before pushing. Do not force-push or overwrite other work.
