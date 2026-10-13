# Private prepared inputs

Prepared sample tables and participant-review records are private and ignored by Git. `samples.tsv` has exact VCF IDs in `sample` and authoritative phenotype labels in `group`.

An optional verified alias table has `metadata_sample` and `vcf_sample` columns. Keep source tables, identity verification, group-assignment audits and provenance together in the private workspace. The preparation script reports unmatched identifiers as unresolved and does not guess aliases. Resolve identity and group assignments before analysis; identity resolution does not establish genotype completeness.

Use a new output path for each preparation or quality-control run. Never publish these participant-level files to the public repository.
