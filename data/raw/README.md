# Private source inputs

Supply a multisample VCF, compatible-build GWAS table and authoritative participant metadata. Configure their paths in `config/settings.json`. Inputs are excluded from the public repository; the private project's release archives contain the supplied data and restoration instructions.

The metadata interface accepts participant identifiers from `name` and source group labels from `Diabetes_Group`. The supplied schema may also include `External Name`. The prepared sample table must contain `sample` and `group` columns and use exact VCF sample IDs. Display `T2D` as Diabetes while preserving the source category in analysis tables.

Use only independently verified one-to-one aliases; never guess participant identity. Preserve original sources and checksums. Verify matching genome builds and allele representations before running the analysis. Credentials must never be placed in these folders or committed.
