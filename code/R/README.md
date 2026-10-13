# Independent R reproduction

These scripts implement the analysis directly in R: authoritative group preparation,
GWAS lead selection, genotype QC, common marker panels, phase-set-aware haplotypes,
group frequencies, participant-level overall/pairwise carrier tests, one panel C,
and an independent audit of saved fine-mapping evidence. They do not run Python.

Use R 4.3 or newer and `jsonlite` (`install.packages("jsonlite")`). For the large
indexed cohort VCF, install `bcftools`; it only extracts the requested windows,
and R performs all filtering and analysis. A small unindexed VCF can be read
directly. PNG/SVG export requires R's Cairo graphics support.

From the repository root, with the local private inputs restored:

```bash
Rscript code/R/test_analysis.R
Rscript code/R/test_audit.R
Rscript code/R/compare_haplotypes.R \
  --genome-build GRCh38 --test-associations \
  --out results/cohort_R_new
Rscript code/R/audit_finemapping.R \
  --run-dir results/cohort_R_new \
  --gwas data/raw/EUR_diagram_t2d_1e-8.tsv \
  --legacy-dir ../GWAS_regions/SuSiE_all_leads \
  --out results/cohort_R_new/variant_evidence
```

The current real-cohort result is **49 Diabetes, 48 Prediabetes, 50 Control**;
19 selected regions fail lead completeness and one has unconfirmed palindromic
orientation. No usable haplotype panels, group association estimates or causal
variants are available. The best pooled lead coverage is 126/147 (85.71%), below
the unchanged 90% requirement. Panel C plots genotype availability; it does not
plot haplotype frequencies. Three image formats are exports of the same figure.

The verified saved run is `results/haplotype_comparison_R_20261013_verified/`.
Its corrected, visually checked figure is in `panel_c_publication/`; its
top-level image exports retain an earlier layout. Current code produces the
corrected layout on new runs. To re-render only the saved aggregate figure:

```bash
Rscript code/R/export_panel_c.R \
  --run-dir results/haplotype_comparison_R_20261013_verified \
  --out results/R_panel_c_new
```

`--test-associations` requests inference but leaves all estimates and p-values
unavailable unless `--reference-genotypes-confirmed` is supported by independent
callability checks. Do not use that flag simply to obtain output. Similarly,
`--gwas-alleles-forward-confirmed` requires actual allele-orientation evidence.
Missing calls are never converted to reference calls. No high-quality 0/0 calls
were observed in the current selected regions, so the provenance gate remains
closed. Unsupported tests retain a reason rather than fabricated estimates.

Settings come from `config/settings.json`; explicit options override them.
Input paths in the settings are relative to the repository; command-line paths
are relative to the working directory. Always choose a new output folder.
To recreate the sample table directly from the final metadata and verified aliases,
add these two arguments to the analysis command:

```text
--metadata data/raw/participant_groups_final.tsv
--id-map data/intermediate/participant_review/verified_id_map.tsv
```

The prepared table is saved locally as `samples_used.tsv`. Neither that file nor
sample-level genotypes/haplotypes should be published publicly. `--samples FILE`
instead uses an already verified exact-VCF-ID sample/group table. Add `--synthetic`
when using software fixtures and `--no-figure` when checking tables only.

## Outputs and methods

The R analysis uses the same main TSV schemas as Python. It writes selected lead
and marker tables, regional QC/allele frequencies, sample haplotypes, group phase
QC, haplotype frequencies, overall comparisons, pairwise comparisons, and
`panel_c_data.tsv`. The sole figure is `panel_c.png` (600 dpi), `panel_c.pdf`, and
`panel_c.svg`. A native `variant_evidence_assessment.tsv` records why causal
resolution is unavailable. The separate audit reads original GWAS and legacy
SuSiE QC/PIP/credible-set/genotype artifacts, writes external lead evidence and
aggregate audit rows, and deliberately writes a header-only candidate table.
This audit only accepts the current unresolved-cohort setting; it does not
pretend to fine-map a future cohort with associated haplotypes.

P-values are compared/ranked using complete decimal digit strings and decimal
exponents, including values smaller than machine precision. Lead clustering uses
consecutive significant-position gaps, not LD. Marker selection is phenotype
blind. Multiple heterozygous sites require phased genotypes and a common explicit
phase set; one heterozygous site or homozygous sites can resolve without it.

Inferential tests use participants as independent units. Overall 2- or 3-group
carrier tables use exact conditional Freeman–Halton enumeration. Pairwise tests
use `stats::fisher.test` p-values, conditional odds ratios and exact 95% intervals
with conditional-likelihood roots refined on log odds to 1e-12;
carrier-frequency differences use Newcombe–Wilson 95% intervals. Overall and
pairwise p-values have separate BH corrections across all loci. All tests remain
exploratory, without adjustment for ancestry, age, sex, batch or relatedness.
Chromosome-copy frequencies are descriptive, not two independent observations
per person. Summary files record settings, input metadata, small-input checksums,
R/package versions, scan scope and interpretation limits.

## Verification

`test_analysis.R` uses R-generated fixtures and tests missing/partial/invalid
genotypes, genotype filtering, multiallelic variants, phase boundaries, phase
flip invariance, p-value precision, exact tests/intervals, BH families, provenance
gates, empty analyses, output preservation and copy/carrier conservation.

For saved Python/R runs, compare full tables and independently audited evidence:

```bash
Rscript code/R/verify_parity.R \
  --python-dir results/haplotype_comparison_20261013/cohort \
  --r-dir results/cohort_R_new \
  --python-audit results/haplotype_comparison_20261013/variant_evidence \
  --r-audit results/cohort_R_new/variant_evidence \
  --out results/cohort_R_new/parity_summary.json
```

The parity summary contains aggregate counts/checks only. Numeric table parity
allows 1e-7 general tolerance because Python formats many outputs to eight
significant figures; conditional odds ratios/CI use the same tolerance after
refining their roots. Exact-test
unit tests use much tighter p-value tolerance. Raw table hashes are reported
honestly: textual floating-point/JSON formatting and gzip headers can differ
even when all analyzed values agree. Indexed R extraction and whole-VCF Python
scanning have different total record counts; all records within selected
windows, their QC, and their downstream results must agree.

The private inputs and legacy artifacts are not part of the public checkout.
