#!/usr/bin/env python3
"""Private descriptive haplotype frequencies and control contrasts, without inference."""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path

from group_labels import display_group

PLOT_FIELDS = ["locus_id", "chrom", "lead_pos", "haplotype_id", "group", "n_group", "n_resolved",
               "haplotype_copies", "resolved_chromosomes", "haplotype_frequency"]
DIFFERENCE_FIELDS = ["locus_id", "chrom", "lead_pos", "haplotype_id", "contrast", "case_group", "control_group",
                     "haplotype_copies_case", "haplotype_copies_control", "resolved_chromosomes_case",
                     "resolved_chromosomes_control", "n_resolved_case", "n_resolved_control", "n_group_case",
                     "n_group_control", "haplotype_frequency_case", "haplotype_frequency_control",
                     "descriptive_copy_frequency_difference", "interpretation_status"]
INTERPRETATION_FIELDS = ["locus_id", "contrast", "haplotype_id_largest_absolute_descriptive_difference",
                         "descriptive_copy_frequency_difference", "n_resolved_case", "n_resolved_control", "interpretation"]
SOURCES = ["run_summary.json", "lead_loci.tsv", "haplotype_group_frequencies.tsv"]
GROUP_ORDER = ["Diabetes", "Prediabetes", "Control"]
CAUTION = ("Descriptive only: frequencies condition on resolved observed chromosomes. Selective missingness and absent "
           "reference calls can distort the contrasts. No association test, confidence interval, adjusted p-value or causal claim is implied.")


def read_table(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_table(path, rows, fields):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({key: "NA" if row[key] is None else row[key] for key in fields} for row in rows)


def load_descriptive(run_dir):
    run = Path(run_dir)
    if not (run / "COMPLETED.txt").is_file():
        raise ValueError("A completed all-region run is required")
    summary = json.loads((run / "run_summary.json").read_text())
    if summary.get("inputs_stable_during_read") is not True:
        raise ValueError("Input stability must be verified")
    if summary.get("haplotype_carrier_tests_performed", 0) or summary.get("haplotype_overall_tests_performed", 0):
        raise ValueError("This descriptive-only export requires a run with no inferential tests")
    leads = read_table(run / "lead_loci.tsv")
    if len(leads) != summary["selected_loci"] or len(leads) != summary["distance_clusters_available"]:
        raise ValueError("The source must include all available GWAS-selected regions")
    if len({r["locus_id"] for r in leads}) != len(leads):
        raise ValueError("Duplicate loci")
    group_map = {}
    for source, n in summary["group_counts"].items():
        label = display_group(source)
        if label in group_map or not isinstance(n, int) or n <= 0:
            raise ValueError("Group counts must be positive and aliases must not be pooled")
        group_map[label] = source
    if set(group_map) != set(GROUP_ORDER):
        raise ValueError("Diabetes, Prediabetes and Control are required")
    group_order = [group_map[g] for g in GROUP_ORDER]
    lead_by_id = {r["locus_id"]: r for r in leads}
    cells, denominators, copy_totals, haps = {}, {}, defaultdict(int), defaultdict(set)
    for source in read_table(run / "haplotype_group_frequencies.tsv"):
        locus, hap, group = (source[k] for k in ["locus_id", "haplotype_id", "group"])
        if locus not in lead_by_id or lead_by_id[locus]["status"] != "analysed" or group not in group_order:
            raise ValueError("Frequency row refers to an ineligible locus or unknown group")
        key = locus, hap, group
        if key in cells:
            raise ValueError("Duplicate haplotype/group frequency row")
        n_group, n_resolved, copies, chromosomes = (int(source[k]) for k in
            ["n_group", "n_resolved", "haplotype_copies", "resolved_chromosomes"])
        if (n_group != summary["group_counts"][group] or not 0 <= n_resolved <= n_group
                or chromosomes != 2 * n_resolved or not 0 <= copies <= chromosomes):
            raise ValueError("Invalid haplotype copy or participant denominator")
        denominator = denominators.setdefault((locus, group), (n_resolved, chromosomes))
        if denominator != (n_resolved, chromosomes):
            raise ValueError("Haplotypes within a group have inconsistent denominators")
        frequency = copies / chromosomes if chromosomes else None
        recorded = source["haplotype_frequency"]
        if frequency is None:
            if recorded not in {"NA", "", "."}:
                raise ValueError("An unobserved group must have an unavailable frequency")
        elif not math.isclose(float(recorded), frequency, abs_tol=1e-7):
            raise ValueError("Haplotype frequency disagrees with chromosome-copy counts")
        lead = lead_by_id[locus]
        cells[key] = dict(locus_id=locus, chrom=lead["chrom"], lead_pos=lead["lead_pos"], haplotype_id=hap,
                          group=group, n_group=n_group, n_resolved=n_resolved, haplotype_copies=copies,
                          resolved_chromosomes=chromosomes, haplotype_frequency=frequency)
        haps[locus].add(hap)
        copy_totals[locus, group] += copies
    if not cells:
        raise ValueError("No descriptive haplotype frequencies are available")
    plotted, differences = [], []
    for lead in leads:
        locus = lead["locus_id"]
        for hap in sorted(haps[locus]):
            for group in group_order:
                if (locus, hap, group) not in cells:
                    raise ValueError("Every local haplotype requires all three group rows")
                plotted.append(cells[locus, hap, group])
            for label in GROUP_ORDER[:2]:
                case, control = group_map[label], group_map["Control"]
                a, b = cells[locus, hap, case], cells[locus, hap, control]
                difference = None if a["haplotype_frequency"] is None or b["haplotype_frequency"] is None else a["haplotype_frequency"] - b["haplotype_frequency"]
                row = dict(locus_id=locus, chrom=lead["chrom"], lead_pos=lead["lead_pos"], haplotype_id=hap,
                           contrast=f"{label} vs Control", case_group=case, control_group=control,
                           descriptive_copy_frequency_difference=difference, interpretation_status="descriptive_only_not_tested")
                for suffix, source in [("case", a), ("control", b)]:
                    row.update({f"{key}_{suffix}": source[key] for key in
                                ["haplotype_copies", "resolved_chromosomes", "n_resolved", "n_group", "haplotype_frequency"]})
                differences.append(row)
        for group in group_order:
            if haps[locus] and copy_totals[locus, group] != denominators[locus, group][1]:
                raise ValueError("Haplotype copies do not sum to the observed chromosomes")
    return summary, plotted, differences


def region_interpretations(differences):
    grouped = defaultdict(list)
    for row in differences:
        grouped[row["locus_id"], row["contrast"]].append(row)
    output = []
    for (locus, contrast), rows in grouped.items():
        available = [r for r in rows if r["descriptive_copy_frequency_difference"] is not None]
        largest = max((abs(r["descriptive_copy_frequency_difference"]) for r in available), default=0)
        tied = [r for r in available if math.isclose(abs(r["descriptive_copy_frequency_difference"]), largest, rel_tol=0, abs_tol=1e-12)]
        best = min(tied, key=lambda r: r["haplotype_id"]) if tied else None
        if best:
            delta = best["descriptive_copy_frequency_difference"]
            text = (f"Largest absolute observed chromosome-copy frequency contrast: {best['haplotype_id']}, {contrast}, "
                    f"{100*delta:+.2f} percentage points among {best['n_resolved_case']} and {best['n_resolved_control']} resolved participants. " + CAUTION)
        else:
            text = "No descriptive contrast is available because a group has no resolved chromosomes. " + CAUTION
        output.append(dict(locus_id=locus, contrast=contrast,
            haplotype_id_largest_absolute_descriptive_difference=best["haplotype_id"] if best else None,
            descriptive_copy_frequency_difference=best["descriptive_copy_frequency_difference"] if best else None,
            n_resolved_case=rows[0]["n_resolved_case"], n_resolved_control=rows[0]["n_resolved_control"], interpretation=text))
    return output


def draw_descriptive(summary, plotted, differences, facets_per_row=3):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.colors import Normalize
    import numpy as np
    if not 1 <= facets_per_row <= 6:
        raise ValueError("facets_per_row must be between 1 and 6")
    loci = list(dict.fromkeys(r["locus_id"] for r in plotted))
    frequencies = {(r["locus_id"], r["haplotype_id"], display_group(r["group"])): r for r in plotted}
    deltas = {(r["locus_id"], r["haplotype_id"], r["contrast"]): r for r in differences}
    haps = {l: list(dict.fromkeys(r["haplotype_id"] for r in plotted if r["locus_id"] == l)) for l in loci}
    cols = min(facets_per_row, len(loci)); nrows = math.ceil(len(loci)/cols)
    heights = [max(2.6, 1.4 + .22 * max(len(haps[l]) for l in loci[i*cols:(i+1)*cols])) for i in range(nrows)]
    figure = plt.figure(figsize=(6*cols, sum(heights)+2.5))
    grid = figure.add_gridspec(nrows, cols, left=.065, right=.975, top=.895, bottom=.060,
                              hspace=.70, wspace=.38, height_ratios=heights)
    frequency_cmap = plt.get_cmap("Blues").with_extremes(bad="#d4d8de")
    difference_cmap = plt.get_cmap("RdBu_r").with_extremes(bad="#d4d8de")
    for index, locus in enumerate(loci):
        ax = figure.add_subplot(grid[index//cols, index%cols])
        patterns = haps[locus]
        matrix = [[frequencies[locus, hap, group]["haplotype_frequency"] for group in GROUP_ORDER] for hap in patterns]
        delta = [[deltas[locus, hap, f"{group} vs Control"]["descriptive_copy_frequency_difference"] for group in GROUP_ORDER[:2]] for hap in patterns]
        data = np.asarray([[np.nan if value is None else value for value in row] for row in matrix])
        contrasts = np.asarray([[np.nan if value is None else value for value in row] for row in delta])
        n = len(patterns)
        ax.imshow(np.ma.masked_invalid(data), cmap=frequency_cmap, norm=Normalize(0, 1), interpolation="nearest", aspect="auto",
                  extent=(-.5, 2.5, n-.5, -.5))
        ax.imshow(np.ma.masked_invalid(contrasts), cmap=difference_cmap, norm=Normalize(-1, 1), interpolation="nearest", aspect="auto",
                  extent=(2.5, 4.5, n-.5, -.5))
        ax.set_xlim(-.5, 4.5); ax.set_ylim(n-.5, -.5)
        for y, values in enumerate([frequency + contrast for frequency, contrast in zip(matrix, delta)]):
            for x, value in enumerate(values):
                label = "NA" if value is None else f"{value:.2f}" if x < 3 else f"{value:+.2f}"
                color = "white" if value is not None and (value > .65 if x < 3 else abs(value) > .65) else "#233445"
                ax.text(x, y, label, ha="center", va="center", fontsize=6.5, color=color)
        labels = []
        for group in GROUP_ORDER:
            row = frequencies[locus, patterns[0], group]
            labels.append(f"{group}\n{row['n_resolved']}/{row['n_group']} people")
        ax.set_xticks(range(5), labels + ["Δ D−C\ncopy freq.", "Δ P−C\ncopy freq."], fontsize=7)
        ax.set_yticks(range(n), [h.removeprefix(locus + "_") for h in patterns], fontsize=7.5)
        ax.tick_params(length=0)
        ax.axvline(2.5, color="#172f46", linewidth=1.1)
        for spine in ax.spines.values():
            spine.set_visible(False)
        first = frequencies[locus, patterns[0], "Diabetes"]
        ax.set_title(f"{locus} · {first['chrom']}:{int(first['lead_pos']):,}\n{n} local patterns", fontsize=10, weight="bold", pad=11)
    scope = "SYNTHETIC FIXTURE · " if summary.get("synthetic") else ""
    figure.text(.055, .985, scope + "Observed haplotype patterns and descriptive contrasts", fontsize=19, weight="bold", color="#17334b", va="top")
    figure.text(.055, .963, f"All {len(loci)} regions with observed patterns · {len(plotted)//3} local haplotypes · No inferential tests", fontsize=11, color="#526579")
    cax1 = figure.add_axes([.10, .932, .30, .008]); cax2 = figure.add_axes([.59, .932, .30, .008])
    colorbar = figure.colorbar(plt.cm.ScalarMappable(norm=Normalize(0, 1), cmap=frequency_cmap), cax=cax1, orientation="horizontal")
    colorbar.set_label("Observed chromosome-copy frequency", fontsize=8); colorbar.ax.tick_params(labelsize=7)
    colorbar = figure.colorbar(plt.cm.ScalarMappable(norm=Normalize(-1, 1), cmap=difference_cmap), cax=cax2, orientation="horizontal")
    colorbar.set_label("Descriptive copy-frequency difference (group − Control)", fontsize=8); colorbar.ax.tick_params(labelsize=7)
    figure.text(.055, .037, "Each frequency uses haplotype copies / (2 × resolved participants). Group labels show resolved / total participants. Gray = unavailable.", fontsize=9, color="#526579")
    figure.text(.055, .023, "Pattern identifiers are local to each region. Selective calls and missing reference genotypes may bias both frequencies and differences.", fontsize=9, color="#526579")
    figure.text(.055, .009, "PRIVATE RESULTS · Descriptive observations only: no evidence of a tested association, progression, protective effect, or causal variant.", fontsize=9, color="#526579")
    return figure


def export_descriptive(run_dir, out, *, facets_per_row=3, dpi=300):
    run, out = Path(run_dir), Path(out)
    hashes = {name: hashlib.sha256((run / name).read_bytes()).hexdigest() for name in SOURCES}
    summary, plotted, differences = load_descriptive(run)
    if out.exists() and any(out.iterdir()):
        raise ValueError("Use a new or empty descriptive output directory")
    if not 1 <= facets_per_row <= 6 or not 72 <= dpi <= 600:
        raise ValueError("Use 1–6 facets per row and 72–600 DPI")
    out.mkdir(parents=True, exist_ok=True)
    write_table(out / "plot_data.tsv", plotted, PLOT_FIELDS)
    write_table(out / "descriptive_frequency_differences.tsv", differences, DIFFERENCE_FIELDS)
    interpretations = region_interpretations(differences)
    write_table(out / "region_interpretations.tsv", interpretations, INTERPRETATION_FIELDS)
    figure = draw_descriptive(summary, plotted, differences, facets_per_row)
    for suffix in ["png", "pdf", "svg"]:
        figure.savefig(out / f"descriptive_haplotypes.{suffix}", dpi=dpi, facecolor="white")
    from matplotlib import pyplot as plt
    plt.close(figure)
    regions = len({r["locus_id"] for r in plotted}); patterns = len(plotted)//3
    text = ["Descriptive haplotype patterns", "", f"{regions} of {summary['selected_loci']} selected regions have observed haplotype patterns; all {patterns} local patterns are shown.",
            "The first three columns show chromosome-copy frequencies among resolved participants. The last two show Diabetes minus Control and Prediabetes minus Control descriptive frequency differences.",
            "No inferential association tests were performed. These displayed differences cannot establish which haplotypes differ in the underlying groups.",
            CAUTION, "Haplotype IDs apply within a locus only; H001 in two regions is not the same chromosome segment.",
            "Largest observed contrasts below are ranked by absolute descriptive differences, not statistical evidence or causal importance.", ""]
    text += [f"{r['locus_id']} — {r['interpretation']}" for r in interpretations]
    text += ["", "No causal variant can be established from these descriptive differences. Validated joint reference/alternate genotypes, phase, independent association analysis and fine-mapping evidence are needed.",
             "All figures, tables and this interpretation are private artifacts."]
    (out / "INTERPRETATION.txt").write_text("\n".join(text) + "\n")
    if hashes != {name: hashlib.sha256((run / name).read_bytes()).hexdigest() for name in SOURCES}:
        raise ValueError("Source analysis changed during export")
    manifest = dict(private_results_only=True, descriptive_only=True, selected_regions=summary["selected_loci"],
                    plotted_regions=regions, local_patterns=patterns, frequency_rows=len(plotted), contrast_rows=len(differences),
                    facets_per_row=facets_per_row, png_dpi=dpi, source_sha256=hashes, synthetic=bool(summary.get("synthetic")))
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--facets-per-row", type=int, default=3)
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()
    print(json.dumps(export_descriptive(args.run_dir, args.out, facets_per_row=args.facets_per_row, dpi=args.dpi), indent=2))
