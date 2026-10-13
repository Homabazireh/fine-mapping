# Run the R code and change options

The R scripts and editable examples are present in both repositories. The public
repository contains code, synthetic tests and methods; real inputs, results and
figures belong in the private project. Restore the required private inputs before
running an analysis and choose a new output directory for each settings change.

## External SNP and indel fine-mapping

Edit [config/external_finemapping.json](../config/external_finemapping.json), then
run from the repository root:

```bash
Rscript code/R/finemap_external.R --config config/external_finemapping.json
```

Set `paths.input_dir` to the prepared directory containing `input_manifest.json`
and the matched variant/LD files. Set `paths.out` to a new results directory.
Relative paths in this configuration resolve from the directory containing the
JSON file. This differs from the older cohort `settings.json`, whose paths resolve
from the project root.

| JSON setting | Command-line override | Default | What changes |
| --- | --- | --- | --- |
| `analysis.variant_types` | `--variant-types both` | `both` | Joint SNP/indel model; use `snps` or `indels` for a restricted model. |
| `analysis.locus` | `--locus LOCUS_ID` | `null` | All prepared regions; specify one identifier for a regional run. |
| `analysis.max_effects` | `--max-effects 5` | `10` | Maximum number of SuSiE single-effect components, capped by variant count. |
| `analysis.coverage` | `--coverage 0.90` | `0.95` | Target posterior mass of model credible sets. |
| `analysis.max_iter` | `--max-iter 300` | `200` | Maximum model iterations; nonconverged results remain unavailable. |
| `analysis.workers` | `--workers 2` | `4` | Number of regional workers. |
| `plot.width`, `plot.height` | `--plot-width 1800 --plot-height 1150` | `1800`, `1150` | Figure dimensions in pixels. |
| `plot.dpi` | `--plot-dpi 150` | `150` | PNG resolution; vector dimensions use pixels divided by DPI. |
| `plot.formats` | `--plot-formats png,pdf,svg` | `["png"]` | Export formats; PDF and SVG support vector editing. |
| `paths.input_dir`, `paths.out` | `--input-dir DIR --out NEW_DIR` | Example paths | Prepared input location and new output location. |

Explicit command-line arguments override configuration values. Command-line paths
resolve from the current working directory. Unknown keys and invalid values are
rejected. Run `Rscript code/R/finemap_external.R --help` for the accepted options.

For example, keep joint SNP/indel fitting while changing model size and credible
set coverage and exporting editable vector figures:

```bash
Rscript code/R/finemap_external.R \
  --config config/external_finemapping.json \
  --max-effects 5 --coverage 0.90 \
  --plot-formats png,pdf,svg \
  --out results/external_settings_comparison
```

Changing `variant_types` changes the model's candidate set. It subsets the
validated variant table and both dimensions of its LD matrix together, preserving
their order. Fewer than two selected variants produces an explicit unavailable
status. Model comparisons use different candidate universes and should not be
interpreted as interchangeable PIPs. `coverage` concerns model credible sets;
the external beta/SE intervals shown in the top panel remain pointwise 95%
normal approximations. Neither setting supplies cohort carrier confidence
intervals or establishes causal variants.

The output records effective settings and configuration provenance alongside the
source and software hashes. Existing runs are preserved. Missing-call,
allele-alignment, reference-provenance and LD validity requirements are not
disabled by these options. See [external methods and preparation](external_finemapping.md)
for the assumptions and required input files.

## RStudio

Open the repository as the working directory, install the required R packages
(`jsonlite`, `openssl`, `susieR`), then run:

```r
source("code/R/finemap_external.R")

run_external_finemapping_config(
  "config/external_finemapping.json",
  overrides = list(
    max_effects = 5,
    coverage = 0.90,
    plot_formats = c("png", "pdf", "svg"),
    out = "results/external_RStudio_new"
  )
)
```

Sourcing the script defines functions without starting an analysis. Override paths
in R are relative to the working directory, as for command-line overrides.

## Cohort analysis and existing pictures

The [R workflow example](../config/R_workflow.example.json) selects an action and
its options. Edit its `action` field and the matching entry in `actions`, then
inspect the command before running it:

```bash
Rscript code/R/run_configured_analysis.R \
  --config config/R_workflow.example.json --dry-run
Rscript code/R/run_configured_analysis.R \
  --config config/R_workflow.example.json
```

Available actions are `cohort_analysis`, `regional_plots`, `haplotype_contrasts`,
`descriptive_haplotypes`, and `external_finemapping`. Only the selected action
runs. Relative workflow paths resolve from the JSON file's directory. A dry run
checks workflow structure and option types and prints the command without creating
results; it can be used before private data are restored. The selected native
script validates model ranges and data contents during an actual run.

For RStudio:

```r
source("code/R/run_configured_analysis.R")
run_configured_analysis("config/R_workflow.example.json", dry_run = TRUE)
run_configured_analysis("config/R_workflow.example.json")
```

For cohort scientific parameters such as flank, marker count and genotype QC,
edit [config/settings.json](../config/settings.json). The wrapper uses that
existing configuration rather than a second copy of its thresholds. Reference
and strand confirmation require independent evidence; they are not casual
workflow switches.

For regional pictures, `include_indels: true` displays SNPs and indels. That
option changes display coverage only. It is distinct from the external model's
`variant_types` setting. The direct R command remains available:

```bash
Rscript code/R/plot_all_regions.R \
  --run-dir results/cohort_R_new \
  --out results/regional_R_new --include-indels
```

See the [R script catalog and tests](../code/R/README.md) for native cohort
analysis, contrast plotting, phase checks and independent Python/R comparison.
