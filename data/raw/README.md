# Raw inputs

Treat this directory as read-only source material. Generate prepared tables
under `../intermediate/` and analysis output under `../../results/`.

- `EUR_diagram_t2d_1e-8.tsv`: default GWAS source; the analysis still applies
  `Pval <= 1e-8` while reading.
- `EUR_diagram_t2d.hg38.tsv`: full GWAS source, available via `--gwas`.
- `participant_groups_final.tsv`: final user-supplied authoritative table with
  columns `External Name`, `Diabetes_Group` and `name`. The grouping column is
  `Diabetes_Group`: 49 T2D, 48 Prediabetes and 50 Control. Participant IDs come
  from `name`. Its provenance is recorded in `participant_groups_final.provenance.json`.
- `participant_categories.tsv`: superseded supplied table, preserved unchanged
  with its `participant_categories.provenance.json` record.
- `participant_metadata.tsv`: earlier supplied table, preserved unchanged.
- `full_LR-WGS.vcf.gz` and `.csi`: relative symbolic links to the originals in
  the parent `Fine_mapping/` directory.
- `group_lists/`: relative symbolic links to the existing `NORM_ID.txt`,
  `Prediabetes_ID.txt` and `T2D_ID.txt` files in `Fine_mapping/`.

The symbolic links avoid copying the cohort files. If this project is moved
away from its parent directory, restore the links or update
`../../config/settings.json` to point to the original inputs. The sequential
VCF reader does not use the `.csi` index.

The final table supersedes the previous table and older group lists. It changes
92 phenotype assignments and changes another 20 labels from `NORM` to `Control`
without changing those participants' displayed phenotype group. All 11
proposed VCF identifier aliases are confirmed; their map is recorded in
`../intermediate/participant_review/verified_id_map.tsv`. Together with 136
exact matches, these account for all 147 participants.
