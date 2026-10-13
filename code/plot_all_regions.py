#!/usr/bin/env python3
"""Export every GWAS-selected region as a separate, private marker figure.

Requires a completed all-cluster run of compare_haplotypes.py (--top-loci 0).
Plots all available PASS biallelic SNPs in each saved lead-centered window,
including windows whose lead/haplotype panel is unusable. No new inference.
"""
from __future__ import annotations

import argparse
from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor
import csv
from datetime import datetime, timezone
import gzip
import hashlib
from html import escape
from itertools import groupby
import json
import math
from pathlib import Path
import re
import shutil
from statistics import median

from comparison_figures import COLORS, PAPER_FIGURE
from group_labels import display_group, ordered_groups

PLOT_FIELDS = ["locus_id", "variant_id", "chrom", "pos", "ref", "alt", "distance_kb",
               "is_exact_lead_match", "group", "n_group", "n_called", "completeness",
               "alt_copies", "called_chromosomes", "observed_alt_fraction"]
INTERPRETATION_FIELDS = ["locus_id", "chrom", "lead_pos", "start", "end", "haplotype_panel_status",
    "group", "n_group", "biallelic_snps", "other_records_omitted", "sites_with_calls", "sites_without_calls",
    "median_completeness", "median_observed_alt_fraction", "lead_n_called", "lead_completeness",
    "lead_observed_alt_fraction", "interpretation"]
CAUTION = ("Observed ALT fractions use only called chromosomes and may be biased by missingness and reference-call "
           "ascertainment. ALT is not automatically a risk allele. These unconnected site measurements do not "
           "show phased haplotypes, test group differences, or identify causal variants.")


def fingerprint(path):
    path = Path(path)
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return {"name": path.name, "bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def read_table(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_run(run_dir, allow_synthetic=False):
    run = Path(run_dir)
    if not (run / "COMPLETED.txt").is_file():
        raise ValueError("Only completed runs can supply regional figures")
    summary = json.loads((run / "run_summary.json").read_text())
    if summary.get("inputs_stable_during_read") is not True:
        raise ValueError("Source inputs were not stable during analysis")
    if summary.get("synthetic") is True and not allow_synthetic:
        raise ValueError("Use --allow-synthetic explicitly for labeled software fixtures")
    leads = read_table(run / "lead_loci.tsv")
    if (not leads or len(leads) != summary.get("selected_loci")
            or len(leads) != summary.get("distance_clusters_available")):
        raise ValueError("All available GWAS clusters are required; rerun the analysis with --top-loci 0")
    if len({r["locus_id"] for r in leads}) != len(leads):
        raise ValueError("Duplicate locus identifiers")
    if not summary.get("group_counts") or any(int(n) <= 0 for n in summary["group_counts"].values()):
        raise ValueError("Positive participant group denominators are required")
    for lead in leads:
        if not re.fullmatch(r"[A-Za-z0-9_-]+", lead["locus_id"]) or not re.fullmatch(r"chr[0-9XYM]+", lead["chrom"]):
            raise ValueError("Unsafe or unsupported locus filename component")
        if not int(lead["start"]) <= int(lead["lead_pos"]) <= int(lead["end"]):
            raise ValueError("Lead lies outside its saved region")
    return summary, leads


def validate_variant(block, lead, totals):
    """Validate both alleles and all groups before exporting count-derived SNP rows."""
    first = block[0]
    definition_keys = ["locus_id", "variant_id", "chrom", "pos", "ref", "alt", "type",
                       "distance_to_lead_bp", "is_exact_lead_match"]
    definition = tuple(first[k] for k in definition_keys)
    pos = int(first["pos"])
    if (first["chrom"] != lead["chrom"] or not int(lead["start"]) <= pos <= int(lead["end"])
            or int(first["distance_to_lead_bp"]) != pos - int(lead["lead_pos"])):
        raise ValueError("Regional position/window/distance mismatch")
    if first["is_exact_lead_match"] not in {"0", "1"} or (
            first["is_exact_lead_match"] == "1" and pos != int(lead["lead_pos"])):
        raise ValueError("Exact lead flag contradicts the saved lead position")
    alleles = [first["ref"], *first["alt"].split(",")]
    cells = {}
    for row in block:
        if tuple(row[k] for k in definition_keys) != definition:
            raise ValueError("Variant definition differs between allele/group rows")
        group = row["group"]
        if group not in totals:
            raise ValueError("Regional row contains an unknown group")
        n, called, chromosomes, copies, index = (int(row[k]) for k in
            ("n_group", "n_called", "called_chromosomes", "allele_copies", "allele_index"))
        if (n != totals[group] or not 0 <= called <= n or chromosomes != 2 * called
                or not 0 <= copies <= chromosomes or not 0 <= index < len(alleles)
                or row["allele"] != alleles[index]):
            raise ValueError("Invalid allele copies or genotype denominators")
        cell = cells.setdefault(group, {"called": called, "chromosomes": chromosomes, "copies": {}})
        if cell["called"] != called or cell["chromosomes"] != chromosomes or index in cell["copies"]:
            raise ValueError("Duplicate allele row or inconsistent genotype denominator")
        cell["copies"][index] = copies
    if set(cells) != set(totals):
        raise ValueError("Variant is missing a participant group")
    for cell in cells.values():
        if set(cell["copies"]) != set(range(len(alleles))) or sum(cell["copies"].values()) != cell["chromosomes"]:
            raise ValueError("Allele copies do not sum to called chromosomes")
    biallelic = (first["type"] == "SNP" and len(alleles) == 2 and all(len(a) == 1 and a in "ACGT" for a in alleles))
    if not biallelic:
        return []
    return [{"locus_id": lead["locus_id"], "variant_id": first["variant_id"], "chrom": first["chrom"],
             "pos": pos, "ref": first["ref"], "alt": first["alt"], "distance_kb": (pos-int(lead["lead_pos"]))/1000,
             "is_exact_lead_match": int(first["is_exact_lead_match"]), "group": group,
             "n_group": totals[group], "n_called": cell["called"], "completeness": cell["called"]/totals[group],
             "alt_copies": cell["copies"][1], "called_chromosomes": cell["chromosomes"],
             "observed_alt_fraction": cell["copies"][1]/cell["chromosomes"] if cell["chromosomes"] else "NA"}
            for group, cell in cells.items()]


class PlotWriters:
    """Bound open descriptors while streaming aggregate rows to every region."""
    def __init__(self, paths, limit=32):
        self.paths, self.limit, self.opened = paths, limit, OrderedDict()

    def write(self, ident, rows):
        if ident not in self.opened:
            if len(self.opened) >= self.limit:
                self.opened.popitem(last=False)[1][0].close()
            path = self.paths[ident]
            fresh = not path.exists()
            handle = gzip.open(path, "at", newline="")
            writer = csv.DictWriter(handle, PLOT_FIELDS, delimiter="\t", lineterminator="\n")
            if fresh:
                writer.writeheader()
            self.opened[ident] = handle, writer
        self.opened.move_to_end(ident)
        self.opened[ident][1].writerows(rows)

    def close(self):
        for handle, writer in self.opened.values():
            handle.close()
        self.opened.clear()


def prepare_plot_data(run_dir, out, summary, leads):
    totals = {g: int(n) for g, n in summary["group_counts"].items()}
    lead_map = {r["locus_id"]: r for r in leads}
    paths, metadata, seen = {}, {}, {r["locus_id"]: set() for r in leads}
    for lead in leads:
        ident = lead["locus_id"]
        folder = out / f"{ident}_{lead['chrom']}_{lead['lead_pos']}"
        folder.mkdir()
        paths[ident] = folder / "plot_data.tsv.gz"
        metadata[ident] = {"lead": lead, "folder": folder.name, "biallelic_snps": 0,
                           "other_records_omitted": 0, "omitted_record_types": {}}
    writers = PlotWriters(paths)
    try:
        with gzip.open(Path(run_dir) / "regional_variant_frequencies.tsv.gz", "rt", newline="") as handle:
            records = csv.DictReader(handle, delimiter="\t")
            for (ident, variant), block in groupby(records, key=lambda r: (r["locus_id"], r["variant_id"])):
                if ident not in lead_map:
                    raise ValueError("Regional row contains an unknown locus")
                if variant in seen[ident]:
                    raise ValueError("Duplicate or non-contiguous regional variant block")
                seen[ident].add(variant)
                block = list(block)
                rows = validate_variant(block, lead_map[ident], totals)
                if rows:
                    metadata[ident]["biallelic_snps"] += 1
                    writers.write(ident, rows)
                else:
                    metadata[ident]["other_records_omitted"] += 1
                    kind = "multiallelic_SNP" if block[0]["type"] == "SNP" and "," in block[0]["alt"] else block[0]["type"]
                    counts = metadata[ident]["omitted_record_types"]
                    counts[kind] = counts.get(kind, 0) + 1
        for ident in paths:
            if not paths[ident].exists():
                writers.write(ident, [])
    finally:
        writers.close()
    return [metadata[r["locus_id"]] for r in leads]


def summarize_region(metadata, rows, summary):
    lead = metadata["lead"]
    records = []
    for group in ordered_groups(summary["group_counts"]):
        cells = [r for r in rows if r["group"] == group]
        calls = [int(r["n_called"]) for r in cells]
        completeness = [float(r["completeness"]) for r in cells]
        alt = [float(r["observed_alt_fraction"]) for r in cells if r["observed_alt_fraction"] != "NA"]
        anchors = [r for r in cells if int(r["is_exact_lead_match"]) == 1]
        anchor = anchors[0] if len(anchors) == 1 else None
        median_call, median_alt = median(completeness) if completeness else "NA", median(alt) if alt else "NA"
        call_text = f"{median_call:.1%}" if median_call != "NA" else "unavailable"
        alt_text = f"{median_alt:.1%}" if median_alt != "NA" else "unavailable"
        text = (f"{display_group(group)}: {len(cells)} SNP sites; {sum(n > 0 for n in calls)} have usable calls. "
                f"Median per-site completeness {call_text}; median ALT-copy fraction among called sites {alt_text}. "
                "These medians describe sites, not independent participant observations. ")
        if anchor:
            text += (f"At the exact SNP lead, {anchor['n_called']}/{summary['group_counts'][group]} participants have usable calls; "
                     f"ALT count is {anchor['alt_copies']}/{anchor['called_chromosomes']} called chromosomes. "
                     + ("Zero called chromosomes means no ALT-fraction estimate. " if int(anchor['called_chromosomes']) == 0 else ""))
        else:
            text += "No unique exact biallelic SNP lead is available for a lead-specific estimate. "
        text += CAUTION
        records.append({"locus_id": lead["locus_id"], "chrom": lead["chrom"], "lead_pos": lead["lead_pos"],
            "start": lead["start"], "end": lead["end"], "haplotype_panel_status": lead["status"], "group": group,
            "n_group": int(summary["group_counts"][group]), "biallelic_snps": metadata["biallelic_snps"],
            "other_records_omitted": metadata["other_records_omitted"], "sites_with_calls": sum(n > 0 for n in calls),
            "sites_without_calls": sum(n == 0 for n in calls), "median_completeness": median_call,
            "median_observed_alt_fraction": median_alt, "lead_n_called": anchor["n_called"] if anchor else "NA",
            "lead_completeness": anchor["completeness"] if anchor else "NA",
            "lead_observed_alt_fraction": anchor["observed_alt_fraction"] if anchor else "NA", "interpretation": text})
    return records


def render_region(task):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    metadata, summary, out, dpi = task
    folder = Path(out) / metadata["folder"]
    rows = read_table(folder / "plot_data.tsv.gz")
    lead = metadata["lead"]
    groups = ordered_groups(summary["group_counts"])
    interpretation = summarize_region(metadata, rows, summary)
    fig, axes = plt.subplots(len(groups), 2, figsize=(13, 3.4 + 1.8*len(groups)), squeeze=False, sharex=True, sharey=True)
    fig.subplots_adjust(left=.105, right=.97, top=.76, bottom=.19, wspace=.17, hspace=.38)
    prefix = "SYNTHETIC FIXTURE | " if summary.get("synthetic") else ""
    fig.text(.055, .956, f"{prefix}{lead['locus_id']} | Regional SNP comparison", fontsize=19, weight="bold", color="#172f46")
    fig.text(.055, .908, f"{lead['chrom']}:{int(lead['start']):,}–{int(lead['end']):,} | "
             f"{summary['genome_build_user_declared']} | GWAS lead {int(lead['lead_pos']):,}\n"
             f"All {metadata['biallelic_snps']:,} available biallelic SNPs shown; "
             f"{metadata['other_records_omitted']:,} other variant records omitted.\n"
             f"Recorded haplotype-panel status: {lead['status']}", fontsize=10, color="#526579", linespacing=1.6, va="top")
    axes[0, 0].set_title("a   Genotype completeness", loc="left", fontsize=12, pad=24)
    axes[0, 1].set_title("b   Observed ALT-copy fraction", loc="left", fontsize=12, pad=24)
    for i, group in enumerate(groups):
        cells = [r for r in rows if r["group"] == group]
        color = COLORS.get(group, "#63758a")
        for j, key in enumerate(["completeness", "observed_alt_fraction"]):
            ax = axes[i, j]
            finite = [r for r in cells if r[key] != "NA"]
            ax.scatter([float(r["distance_kb"]) for r in finite], [100*float(r[key]) for r in finite],
                       s=8, alpha=.55, color=color, linewidths=0, rasterized=True)
            anchors = [r for r in cells if int(r["is_exact_lead_match"]) == 1]
            if len(anchors) == 1:
                numerator, denominator = ((anchors[0]["n_called"], summary["group_counts"][group]) if j == 0 else
                                          (anchors[0]["alt_copies"], anchors[0]["called_chromosomes"]))
                label = f"Lead: {numerator}/{denominator}" if int(denominator) else "Lead ALT fraction: NA (zero calls)"
            else:
                label = "Lead: no unique exact SNP estimate"
            ax.text(1, 1.04, label, transform=ax.transAxes, ha="right", fontsize=8, color="#526579")
            if len(anchors) == 1 and anchors[0][key] != "NA":
                ax.scatter([float(anchors[0]["distance_kb"])], [100*float(anchors[0][key])], marker="D", s=54,
                           color=color, edgecolor="#172f46", linewidth=.8, zorder=4)
            ax.axvline(0, color="#a53d33", ls="--", lw=.8)
            ax.set_xlim((int(lead["start"])-int(lead["lead_pos"]))/1000, (int(lead["end"])-int(lead["lead_pos"]))/1000)
            ax.set_ylim(-4, 104)
            ax.set_yticks([0, 25, 50, 75, 100])
            ax.grid(axis="y", color="#e8edf2", lw=.7)
            ax.spines[["top", "right"]].set_visible(False)
            if not finite:
                ax.text(.5, .5, "No SNP sites" if not cells else "No called chromosomes; ALT fraction unavailable",
                        transform=ax.transAxes, ha="center", va="center", fontsize=9, color="#526579")
            if j == 1:
                ax.text(.01, .04, f"{sum(int(r['n_called']) == 0 for r in cells)} sites without calls: ALT fraction NA",
                        transform=ax.transAxes, fontsize=8, color="#526579",
                        bbox={"facecolor": "white", "edgecolor": "none", "alpha": .85})
        axes[i, 0].set_ylabel(f"{display_group(group)}\nn = {summary['group_counts'][group]}\nPercent", fontsize=10, color=color)
    for ax in axes[-1]:
        ax.set_xlabel("Distance from GWAS lead (kb)", fontsize=10, labelpad=9)
    fig.text(.055, .125, "Left: usable participants / all participants in each group. Right: ALT copies / (2 × usable participants).\n"
             "Every plotted site is a SNP; sites with zero calls have 0% completeness and no ALT-fraction estimate.\n"
             "Missingness and reference-call ascertainment can bias observed fractions. No phase, tested group difference or causal inference.\n"
             "Dashed line: GWAS lead position. Outlined diamond: unique exact SNP lead when available. ALT is not necessarily a risk allele.",
             fontsize=8.5, color="#526579", va="top", linespacing=1.65)
    for suffix in ("png", "pdf"):
        fig.savefig(folder / f"region.{suffix}", dpi=dpi, facecolor="white")
    plt.close(fig)
    panel_context = ("This region yielded a descriptive haplotype panel in the saved analysis. Its separate haplotype frequencies "
                     "are conditional on callable, phase-resolved participants and do not establish a group association."
                     if lead["status"] == "analysed" else
                     "This region did not yield an eligible haplotype panel under the saved genotype and orientation requirements.")
    note = [f"# {lead['locus_id']}: regional marker interpretation", "",
        f"Window: {lead['chrom']}:{lead['start']}–{lead['end']} ({summary['genome_build_user_declared']}).",
        f"All {metadata['biallelic_snps']} available biallelic SNPs are included; {metadata['other_records_omitted']} other records are omitted.",
        "Omitted record classes: " + (", ".join(f"{k}: {v}" for k, v in metadata["omitted_record_types"].items()) or "none") + ".",
        f"Haplotype-panel status from the unchanged analysis: `{lead['status']}`.", "",
        panel_context, "",
        "![Regional marker comparison](region.png)", "", "[Vector PDF](region.pdf) · [Plotted counts](plot_data.tsv.gz)", "",
        *[r["interpretation"] + "\n" for r in interpretation],
        "A missing or unusable lead never suppresses neighboring SNPs. Markers are individual sites; point spacing represents genomic distance.",
        "Group medians are unweighted summaries across sites; they are not tests or population allele-frequency estimates.",
        f"Layout inspiration: [Figueroa et al., Figure 1]({PAPER_FIGURE}).", ""]
    (folder / "INTERPRETATION.md").write_text("\n".join(note))
    files = {name: fingerprint(folder / name) for name in ["region.png", "region.pdf", "plot_data.tsv.gz", "INTERPRETATION.md"]}
    return {"locus_id": lead["locus_id"], "folder": metadata["folder"], "biallelic_snps": metadata["biallelic_snps"],
            "other_records_omitted": metadata["other_records_omitted"], "omitted_record_types": metadata["omitted_record_types"],
            "files": files, "interpretations": interpretation}


def export_all_regions(run_dir, out, *, dpi=300, workers=4, allow_synthetic=False):
    run, out = Path(run_dir), Path(out)
    summary, leads = load_run(run, allow_synthetic)
    if out.exists() and any(out.iterdir()):
        raise ValueError("Use a new or empty regional output directory")
    if not 72 <= dpi <= 600 or workers < 1:
        raise ValueError("DPI must be 72–600 and workers positive")
    source_paths = [run / name for name in ["run_summary.json", "lead_loci.tsv", "regional_variant_frequencies.tsv.gz"]]
    source_before = {path.name: fingerprint(path) for path in source_paths}
    out.mkdir(parents=True, exist_ok=True)
    metadata = prepare_plot_data(run, out, summary, leads)
    tasks = [(item, summary, str(out), dpi) for item in metadata]
    results = []
    if workers == 1:
        iterator = map(render_region, tasks)
        executor = None
    else:
        executor = ProcessPoolExecutor(max_workers=workers)
        iterator = executor.map(render_region, tasks)
    try:
        for i, result in enumerate(iterator, 1):
            results.append(result)
            if i % 25 == 0 or i == len(tasks):
                print(f"Rendered {i}/{len(tasks)} regions", flush=True)
    finally:
        if executor:
            executor.shutdown()
    if source_before != {path.name: fingerprint(path) for path in source_paths}:
        raise ValueError("Source analysis changed during export; discard incomplete output")
    with (out / "interpretations.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, INTERPRETATION_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for result in results:
            writer.writerows(result.pop("interpretations"))
    shutil.copyfile(run / "regional_variant_frequencies.tsv.gz", out / "all_regional_variant_frequencies.tsv.gz")
    intro = (f"All {len(leads)} GWAS-selected lead-centered windows are included, not only the previous top loci. "
             "Each region has its own picture and explanation. All available PASS biallelic SNPs are plotted, "
             "including regions without usable haplotype panels. Multiallelic SNPs and indels are counted as omissions "
             "and remain available in the complete regional allele-count table. "
             f"The source analysis produced {summary.get('locus_status_counts', {}).get('analysed', 0)} descriptive haplotype panels; "
             "these marker figures do not substitute for the separate haplotype-frequency tables. " + CAUTION)
    md = ["# All-region marker comparisons", "", intro, "", "[Region/group interpretation table](interpretations.tsv)", "",
          "[Complete regional allele-count table, including omitted variant types](all_regional_variant_frequencies.tsv.gz)", "",
          "| Region | Location | SNPs | Interpretation | PNG | PDF |", "| --- | --- | ---: | --- | --- | --- |"]
    html_rows = []
    for result, lead in zip(results, leads):
        folder = result["folder"]
        md.append(f"| {lead['locus_id']} | {lead['chrom']}:{lead['start']}–{lead['end']} | {result['biallelic_snps']} | "
                  f"[Read]({folder}/INTERPRETATION.md) | [PNG]({folder}/region.png) | [PDF]({folder}/region.pdf) |")
        html_rows.append(f"<tr><td>{escape(lead['locus_id'])}</td><td>{escape(lead['chrom'])}:{lead['start']}–{lead['end']}</td>"
            f"<td>{result['biallelic_snps']}</td><td><a href='{folder}/region.png'>PNG</a> · <a href='{folder}/region.pdf'>PDF</a> · "
            f"<a href='{folder}/INTERPRETATION.md'>Interpretation</a></td></tr>")
    (out / "index.md").write_text("\n".join(md) + "\n")
    (out / "index.html").write_text("<!doctype html><meta charset='utf-8'><title>All-region marker comparisons</title>"
        "<style>body{font:16px system-ui;max-width:1100px;margin:40px auto;padding:0 20px;color:#213547}"
        "table{border-collapse:collapse;width:100%}td,th{padding:8px;border-bottom:1px solid #ddd;text-align:left}</style>"
        "<h1>All-region marker comparisons</h1><p>" + escape(intro) + "</p><p><a href='interpretations.tsv'>Interpretation table</a> · "
        "<a href='all_regional_variant_frequencies.tsv.gz'>Complete regional allele-count table</a></p>"
        "<table><tr><th>Region</th><th>Window</th><th>SNPs</th><th>Files</th></tr>" + "".join(html_rows) + "</table>")
    manifest = {"schema_version": 1, "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "private_results_only": True, "synthetic": bool(summary.get("synthetic", False)), "source_run": run.name,
        "selection": "all_available_GWAS_distance_clusters", "n_regions": len(leads),
        "n_figures": len(results), "png_dpi": dpi, "group_counts": summary["group_counts"],
        "haplotype_panel_status_counts": summary.get("locus_status_counts", {}),
        "genome_build": summary["genome_build_user_declared"], "source_sha256": source_before,
        "parameters": {k: summary["parameters"].get(k) for k in
                       ["p_threshold", "lead_gap", "top_loci", "flank", "min_dp", "min_gq", "min_call_rate"]},
        "interpretation": CAUTION, "regions": results}
    manifest["index_files"] = {name: fingerprint(out / name) for name in
                               ["interpretations.tsv", "index.md", "index.html", "all_regional_variant_frequencies.tsv.gz"]}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    (out / "COMPLETED.txt").write_text(f"Completed {len(results)} separate regional marker figures. Results and interpretations are private only.\n")
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--allow-synthetic", action="store_true")
    args = parser.parse_args()
    export_all_regions(args.run_dir, args.out, dpi=args.dpi, workers=args.workers, allow_synthetic=args.allow_synthetic)


if __name__ == "__main__":
    main()
