# Raw inputs

Treat this directory as read-only source material. Generate prepared tables
under `../intermediate/` and analysis output under `../../results/`.

- `EUR_diagram_t2d_1e-8.tsv`: default GWAS source; the analysis still applies
  `Pval <= 1e-8` while reading.
- `EUR_diagram_t2d.hg38.tsv`: full GWAS source, available via `--gwas`.
- `participant_categories.tsv`: latest user-supplied authoritative table; exact
  categories T2D, Prediabetes and NORM. Its provenance is recorded in
  `participant_categories.provenance.json`.
- `participant_metadata.tsv`: earlier supplied table, preserved unchanged.
- `full_LR-WGS.vcf.gz` and `.csi`: relative symbolic links to the originals in
  the parent `Fine_mapping/` directory.
- `group_lists/`: relative symbolic links to the existing `NORM_ID.txt`,
  `Prediabetes_ID.txt` and `T2D_ID.txt` files in `Fine_mapping/`.

The symbolic links avoid copying the cohort files. If this project is moved
away from its parent directory, restore the links or update
`../../config/settings.json` to point to the original inputs. The sequential
VCF reader does not use the `.csi` index.

The authoritative category table supersedes the older group lists. Eleven
VCF identifier aliases still require confirmation; see
`../intermediate/participant_review/proposed_id_map.tsv`.
