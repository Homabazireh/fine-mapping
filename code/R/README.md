# Independent R reproduction

**[Editable R options, command-line examples and RStudio entry points](../../docs/R_options.md)** · [External model settings](../../config/external_finemapping.json) · [Workflow action settings](../../config/R_workflow.example.json).

These scripts implement the analysis directly in R: authoritative group preparation,
GWAS lead selection, genotype QC, common marker panels, phase-set-aware haplotypes,
group frequencies, participant-level overall/pairwise carrier tests, panel C,
and an independent audit of saved fine-mapping evidence. They do not run Python.

Use R 4.3 or newer and `jsonlite` (`install.packages("jsonlite")`). For the large
indexed cohort VCF, install `bcftools`; it only extracts the requested windows,
and R performs all filtering and analysis. A small unindexed VCF can be read
directly. PNG/SVG export requires R's Cairo graphics support.

From the repository root, with the local private inputs restored:

```bash
Rscript code/R/test_analysis.R
Rscript code/R/test_audit.R
Rscript code/R/test_haplotype_contrasts.R
Rscript code/R/test_all_regions.R
Rscript code/R/test_descriptive_haplotypes.R
Rscript code/R/test_parallel_analysis.R
Rscript code/R/test_external_finemapping.R
Rscript code/R/test_configured_analysis.R
Rscript code/R/compare_haplotypes.R \
  --genome-build GRCh38 --test-associations \
  --out results/cohort_R_new
Rscript code/R/audit_finemapping.R \
  --run-dir results/cohort_R_new \
  --gwas data/raw/EUR_diagram_t2d_1e-8.tsv \
  --legacy-dir ../GWAS_regions/SuSiE_all_leads \
  --out results/cohort_R_new/variant_evidence
```

Cohort figures, tables, counts and interpretations are private. Restore them from
the private [fine-mapping-2 project](https://github.com/Homabazireh/fine-mapping-2).
Panel C plots genotype availability, with PNG/PDF/SVG exports of the same figure.
To re-render only the saved aggregate figure:

```bash
Rscript code/R/export_panel_c.R \
  --run-dir results/cohort_R_new \
  --out results/R_panel_c_new
```

`--test-associations` requests inference but leaves all estimates and p-values
unavailable unless `--reference-genotypes-confirmed` is supported by independent
callability checks. Do not use that flag simply to obtain output. Similarly,
`--gwas-alleles-forward-confirmed` requires actual allele-orientation evidence.
Missing calls are never converted to reference calls. Unsupported tests retain
the relevant quality or provenance reason rather than fabricated estimates.

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
`panel_c_data.tsv`. Panel C is exported as `panel_c.png` (600 dpi), `panel_c.pdf`, and
`panel_c.svg`. A native `variant_evidence_assessment.tsv` records why causal
resolution is unavailable. The separate audit reads original GWAS and legacy
SuSiE QC/PIP/credible-set/genotype artifacts, writes external lead evidence and
aggregate audit rows, and deliberately writes a header-only candidate table.
This audit only accepts an unresolved cohort setting; it does not
pretend to fine-map a future cohort with associated haplotypes.

The requested disease-control figure is generated separately from saved tables:

```bash
Rscript code/R/plot_haplotype_contrasts.R \
  --run-dir results/cohort_R_new --out results/cohort_R_new/contrasts
```

It compares Diabetes versus Control and Prediabetes versus Control across every
available panel. The overview shows observed participant-carrier differences;
`index.md` links separate regional pictures with carrier/resolved denominators
and interpretations. Descriptive-only differences are computed from counts,
retain their blocked `inference_status`, and have no confidence intervals or
p/q values. Only valid tested rows display intervals and the original full-family
BH q-values. Reversed comparisons are reoriented. Unavailable observations
remain missing, and actual observed zero differences remain zero.
`contrast_data.tsv` keeps descriptive carrier differences separate from
inferential fields and chromosome-copy frequencies. `INTERPRETATION.txt`
explains the effects of missingness and the limits of causal claims.
All figures, tables and interpretations stay in the private project.

For separate pictures of every GWAS region, first select every available cluster
with `--top-loci 0`, then run the independent regional plotter on the completed
aggregate run:

```bash
Rscript code/R/compare_haplotypes.R --genome-build GRCh38 \
  --top-loci 0 --test-associations --no-figure --out results/all_regions_R_new
Rscript code/R/plot_all_regions.R \
  --run-dir results/all_regions_R_new --out results/all_regions_R_figures_new \
  --include-indels
```

For a large indexed VCF on Unix, use the parallel native R wrapper instead of
the sequential analysis command. Choose a worker count suited to available CPU
and memory; overlapping windows stay together and the final association/BH step
runs once across all regions:

```bash
Rscript code/R/compare_haplotypes_parallel.R --workers 8 \
  --genome-build GRCh38 --top-loci 0 --test-associations --no-figure \
  --out results/all_regions_R_parallel_new
```

This wrapper uses the same native R QC/phasing/statistical code and requires
`bcftools` plus an indexed VCF. `test_parallel_analysis.R` checks overlapping
windows, chromosome boundaries, aggregate tables and scan counts against the
sequential implementation on synthetic data.

The plotter also accepts the identical aggregate tables from a completed Python
all-cluster run, avoiding a repeated raw scan when only independent plotting is
needed. The manifest identifies the source implementation. It refuses a run
that includes only a subset of available clusters. Each region gets its own
PNG/PDF, complete plotting table and interpretation, including regions with
failed leads or no biallelic SNPs. A searchable browser index links every region.

For graphics-only changes, regenerate pictures in parallel from the already
validated saved plot tables. This checks source hashes and copies all numerical
tables and interpretations unchanged into a new directory:

```bash
Rscript code/R/replot_all_regions.R \
  --run-dir results/all_regions_R_new \
  --plot-dir results/all_regions_R_figures_new \
  --out results/all_regions_R_figures_redrawn --workers 8
```

With `--include-indels`, the panels show available PASS biallelic SNPs and
sequence-resolved small indels by group, using circles and triangles respectively.
Separate class counts and `variant_type` are retained in output tables/manifests.
Without this option, the earlier SNP-only schema remains available.
Genotype completeness
uses all participants in that group, and observed ALT fractions use twice the
number of called participants. REF plus ALT copies are independently checked
against that denominator. A site with zero calls has unavailable ALT frequency,
while its completeness is zero. Medians respect those different denominators.
Multiallelic and unsupported records have explicit omission counts and remain in the
full aggregate variant table. Allele-copy counts describe diploid genotypes,
not copy-number variation. These descriptive pictures do not identify
haplotype associations or causal variants. Results and interpretations are
private only.

Eligible panels can still have descriptive observed haplotype frequencies when
association tests are unavailable. Plot every observed local pattern, with
three group frequencies and two descriptive disease-minus-control columns:

```bash
Rscript code/R/plot_descriptive_haplotypes.R \
  --run-dir results/all_regions_R_new --out results/descriptive_haplotypes_R_new
```

The plot uses chromosome-copy frequencies with denominator twice the resolved
participant count, and labels each region's phase-resolution denominators.
All local haplotype patterns are retained. Gray cells mean unavailable, and
missing frequencies remain missing in both contrast columns. Separate common
color scales distinguish frequencies from descriptive differences. No
confidence intervals, p-values, q-values or significance marks are attached to
these copy-frequency contrasts. They are conditional on observed, callable,
phase-resolved participants and can be biased by missingness. The accompanying
tables report every pattern and the largest absolute descriptive contrast per
region; that ranking does not identify an associated or causal haplotype.
These outputs and interpretations are private.

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
  --python-dir results/cohort_Python_new \
  --r-dir results/cohort_R_new \
  --python-audit results/cohort_Python_new/variant_evidence \
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

## External GWAS fine-mapping

The separate `finemap_external.R` workflow requires `susieR`, `jsonlite` and
`openssl`. It consumes checksummed, ALT-oriented summary statistics and signed
LD prepared with NumPy and `bcftools`. See the [complete extraction, reference
acquisition and fitting instructions](../../docs/external_finemapping.md).
Use `Rscript code/R/test_external_finemapping.R` for synthetic regression checks.
Its posterior probabilities and credible sets describe an exploratory external
study model; they do not establish cohort carrier associations or validated
causal variants.
