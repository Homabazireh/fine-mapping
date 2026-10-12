# Prepared data and review tables

Store derived input tables here. Raw sources belong in `../raw/`.

`participant_review/summary.json` records your authoritative category source.
The 88 historical differences with the old ID lists are resolved in favour of
your table. Eleven proposed VCF aliases remain pending in `proposed_id_map.tsv`.
The earlier comparison files are retained as history.

`cohort_preflight/` records an indexed genotype-completeness check at the
selected GWAS lead sites, independent of participant category matching.

The configured final analysis input is `samples.tsv`, with `sample` and `group`
columns and verified VCF IDs. While the aliases are pending,
`samples_from_authoritative_categories_for_review.tsv` contains 136 exact
assignments and 11 UNRESOLVED rows. No final sample table has been created yet.

`demo/` is created when running the synthetic demonstration. Its participant
labels are generated examples and do not apply to the cohort.
