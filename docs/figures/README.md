# Panel c: usable genotypes across 20 regions

![Panel c: usable lead genotypes in Diabetes, Prediabetes and Control, with pooled counts and a 90 percent threshold](panel_c.png)

[PNG: 600 dpi](panel_c.png) · [PDF: scalable](panel_c.pdf) · [SVG: editable](panel_c.svg)

[Read the plain-language report](REPORT.txt) · [Download the self-contained HTML report](report.html)

The PNG is **7,800 × 6,240 pixels**. PDF and SVG preserve sharp text and markers
when enlarged. Only panel c is included in this figure page and its report.

## What the figure means

The plot uses **147 real participants: 49 Diabetes, 48 Prediabetes and 50
Control**. Each row is a selected DNA region's GWAS lead variant, not a person.
The labels L0001–L0020 follow the public GWAS ranking; chromosome and GRCh38
position are shown alongside them.

Colored symbols show the percentage of each group with a usable genotype.
Black diamonds combine all participants. The counts on the right show usable
calls out of 147. The dashed line is the **90% pooled completeness requirement**.
A point further right means more complete genotype data, not greater diabetes
risk. The highlighted L0010 has the best completeness: **126/147, or 85.7%**.

**None of the 20 leads meets 90%.** The completed run produced zero usable
haplotype panels and zero haplotype association tests. The figure therefore
cannot establish differences or equivalence between the groups' haplotypes.
The separate preliminary VCF check was skipped; the main analysis retained
its recorded genotype and phase requirements. Drawing this figure only reads
saved result tables and does not rescan the VCF.

The layout is inspired by [Figueroa et al. (2024), Figure 1](https://www.nature.com/articles/s41588-024-01719-5/figures/1).
All plotted values come from this cohort. This is a genotype-completeness
figure, not a replication of the paper's founder-haplotype analysis.

## Files and reproduction

[Plotted aggregate counts](panel_c_data.tsv) · [Provenance and checksums](panel_c_manifest.json)

The full local output is saved at `results/cohort_no_preflight/panel_c/`.
To regenerate it into a new folder from the project root:

```bash
.venv/bin/python code/export_panel_c.py \
  --run-dir results/cohort_no_preflight \
  --out results/cohort_no_preflight/panel_c_new
```

This command requires the completed local run's tables. A GitHub clone includes
the published figure and report, but does not include the raw genomic data,
participant tables or full local run. See the [analysis guide](../analysis.md)
and [literature comparison](../../literature_review/figueroa_2024_zfhx3_sca4/review.md).
