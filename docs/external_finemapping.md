# External SNP and indel fine-mapping

This workflow fits exploratory models from external GWAS summary statistics and
public reference genotypes. Cohort carrier associations remain a separate
analysis with their existing genotype and phase requirements. Preserve a new
output directory at every stage; keep region selections, results, figures and
interpretations in the private project.

## Change R model and figure options

Use the [editable options guide](R_options.md) and
[config/external_finemapping.json](../config/external_finemapping.json) to select
variant classes, model size, credible-set coverage, iterations, workers and
PNG/PDF/SVG exports. Configuration is optional; existing command-line calls
retain their scientific defaults. Every run records its effective settings and
requires a new output directory.

```bash
Rscript code/R/finemap_external.R --config config/external_finemapping.json
```

## Inputs and dependencies

Use Python 3.12 or newer (tested with Python 3.14). Install the Python packages
with `python -m pip install -r requirements.txt`, including NumPy for
signed LD calculations. The cohort VCF reader and external statistics extractor
use the Python standard library; that statement does not cover LD preparation
or plotting. Install `bcftools`, an indexed GRCh38 reference FASTA, and the R
packages `susieR`, `jsonlite` and `openssl`. `samtools faidx` can create the FASTA
index if needed.

The extractor currently supports the explicitly selected
`ukb-topmed-t2d-250.2` profile. Place a verified local copy of the full
[PheWeb T2D summary-statistic download](https://pheweb.org/UKB-TOPMed/download/250.2)
at the input path. The portal documents GRCh38 coordinates and ALT-oriented
allele frequencies/effects in its [source description](https://pheweb.org/UKB-TOPMed/about);
trait-level case/control totals come from the [T2D phenotype page](https://pheweb.org/UKB-TOPMed/pheno/250.2).
The profile declaration and local checksum do not prove byte identity to a
remote release. Preserve the download record and verify that the supplied file
matches this source, trait and build before using the profile.

The window table needs unique, filename-safe `locus_id` values and columns
`chrom`, `lead_pos`, `start`, `end`. Use `chr1` through `chr22` and positive
one-based inclusive coordinates with each lead inside its window. The extractor
and reference fetcher currently reject overlapping windows explicitly; they do
not silently assign a shared variant to one of several regions. Saved
lead-centered windows need not cover a complete fine-mapping locus. Review
boundaries and run expanded, nonoverlapping loci when the scientific question
requires them.

## Extract and acquire reference genotypes

Run commands from the repository root, adjusting these generic input paths:

```bash
.venv/bin/python code/extract_external_gwas.py \
  --gwas data/raw/external_t2d.tsv.gz \
  --leads results/all_regions_run/lead_loci.tsv \
  --source-profile ukb-topmed-t2d-250.2 \
  --out results/external_gwas_extracted
.venv/bin/python code/fetch_external_ld_reference.py \
  --leads results/external_gwas_extracted/lead_loci.tsv \
  --fasta data/raw/GRCh38.fa \
  --out results/external_ld_reference --workers 4
```

Extraction accepts plain or gzip-compressed tab-separated input. It retains
eligible SNPs and unequal-length ACGT allele sequences with stored allele
lengths up to `--max-indel-length + 1`, preserving source strings and precision.
It applies no association p-value threshold. MNVs, symbolic/multiallelic alleles,
longer indels and invalid statistics are counted as exclusions. Duplicate exact
variants inside the selected windows stop extraction, even when one duplicate
has invalid statistics. Empty regions retain header-only files. The exact lead
table is copied automatically for the next stage. No completion marker is
written after a failed scan.

The reference fetcher uses the named
[IGSR high-coverage phased release](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/)
and its documented original unrelated European sample subset. It queries the
requested windows, selects biallelic SNP/indel records, and normalizes against
the supplied FASTA with fatal reference mismatches. Acceptance of retained
`FILTER=.` is confined to this already-filtered public release; it does not
change cohort PASS requirements. Source URLs, sample order, commands, hashes,
and before/after identity checks for the window table, FASTA and index are
recorded. Reference acquisition is not cohort genotype recovery.

## Prepare aligned LD and fit

```bash
.venv/bin/python code/prepare_external_finemapping.py \
  --source-dir results/external_gwas_extracted \
  --reference-dir results/external_ld_reference \
  --fasta data/raw/GRCh38.fa \
  --out results/external_finemapping_prepared --workers 4
Rscript code/R/finemap_external.R \
  --input-dir results/external_finemapping_prepared \
  --out results/external_finemapping_fit --workers 4
```

Both preparation and fitting support `--locus LOCUS_ID` for an explicitly
limited run. Omitting that option processes every available input region.
Preparation records each requested region and all observable exclusion
reasons. A completion marker indicates processing finished; inspect statuses
to determine which regions have usable model inputs or fitted results.

Preparation normalizes external GWAS alleles, requires exact normalized
REF/ALT matches, and excludes ambiguous duplicate keys. Reference genotypes
must be completely observed and polymorphic: missing calls are not imputed.
Signed correlations use ALT dosages in an explicit shared variant order.
Common-variant MAF and GWAS/reference ALT-frequency checks further restrict the
available variants. These exclusions limit what the fitted model can explain.

The R model validates hashes, variant order and correlation-matrix properties,
then uses SuSiE-RSS with ALT beta/SE statistics and external LD. It records LD/z
diagnostics, convergence, model settings, SNP/indel classes and window-boundary
flags. The primary binary-trait approximation uses balanced case/control
effective sample size, with nominal total sample size as a sensitivity fit.
Neither is a measured per-variant sample size. Diagnostic thresholds are
reporting heuristics; severe inconsistencies or failed fits leave posterior
probabilities unavailable. Review both the primary and sensitivity outputs.

See the [SuSiE-RSS reference](https://stephenslab.github.io/susieR/reference/susie_rss.html)
and [LD/summary-statistic diagnostics](https://stephenslab.github.io/susieR/articles/susierss_diagnostic.html)
for the underlying model and diagnostic context. Posterior inclusion
probabilities and credible sets remain conditional on the retained variants,
LD, model and window. They are not validated causal-variant sets, and they do
not establish Diabetes/Prediabetes/Control differences in the cohort.

## Checks

```bash
.venv/bin/python code/run_tests.py
Rscript code/R/test_external_finemapping.R
git diff --check
```

Extractor and reference-acquisition tests use small synthetic files and mocked
network/process boundaries; they do not download a reference panel. LD
preparation tests include synthetic `bcftools` normalization and missingness
checks. R tests validate model inputs and reporting behavior on synthetic data.
Keep private source manifests and checksums with each completed run.
