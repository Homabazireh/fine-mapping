# Prepared data and review tables

Store derived input tables here. Raw sources belong in `../raw/`.

`participant_review/summary.json` records the authoritative
`data/raw/participant_groups_final.tsv` source: 49 T2D, 48 Prediabetes and 50
Control. Compared with the previous table, 92 participants change phenotype
group and another 20 change only the category label from `NORM` to `Control`.
All 11 proposed VCF aliases are confirmed in
`participant_review/verified_id_map.tsv`. Together with 136 exact matches,
these account for all 147 participants. The earlier review is retained in
the archive recorded by the summary.

`cohort_preflight/` records an indexed genotype-completeness check at the
selected GWAS lead sites, independent of participant category matching.
`confirmed_groups_qc/summary.json` and `group_lead_genotype_qc.tsv` record the
current check using the final table and verified map: 49 Diabetes, 48 Prediabetes
and 50 Control participants. The exact counts are 45/44/47 and verified aliases
are 4/4/3, respectively. Select this summary explicitly with
`build_report.py --preflight`. The earlier `final_groups_qc/` check used only
the 136 exact matches and is retained as history.

The configured final analysis input is `samples.tsv`, with `sample` and `group`
columns and verified VCF IDs. All 147 participants now have assignments in this
final table; none has an unresolved identifier. The earlier
`samples_from_authoritative_categories_for_review.tsv` candidate is retained
for source review and is not the configured final analysis input.

`demo/` is created when running the synthetic demonstration. Its participant
labels are generated examples and do not apply to the cohort.
