#!/usr/bin/env python3
"""Render private Diabetes/Prediabetes versus Control haplotype comparisons.

Unavailable comparisons are status cells, never fabricated zero effects.
The figure consumes aggregate completed-run tables; it does not rerun inference.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path

from group_labels import display_group


FIELDS = ["locus_id", "chrom", "lead_pos", "haplotype_id", "contrast", "case_group", "control_group",
          "panel_status", "comparison_status", "reason", "n_resolved_case", "n_resolved_control",
          "carriers_case", "carriers_control", "carrier_frequency_case", "carrier_frequency_control",
          "carrier_frequency_difference", "CI95_low", "CI95_high", "p_fisher_two_sided",
          "q_BH_all_pairwise_tests", "haplotype_frequency_case_descriptive", "haplotype_frequency_control_descriptive"]
SOURCES = ["run_summary.json", "lead_loci.tsv", "haplotype_comparisons.tsv", "haplotype_group_frequencies.tsv"]


def read_table(path):
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(f"Missing or duplicate table headers: {path}")
        return list(reader)


def numeric(value, field, *, missing=False):
    if value in {None, "", "NA", "."}:
        if missing:
            return None
        raise ValueError(f"Tested comparison is missing {field}")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"Nonfinite {field}")
    return result


def load_contrasts(run_dir):
    """Normalize both requested contrasts and retain the original BH family."""
    run_dir = Path(run_dir)
    if not (run_dir / "COMPLETED.txt").is_file():
        raise ValueError("Only completed analysis runs can be plotted")
    summary = json.loads((run_dir / "run_summary.json").read_text())
    if summary.get("inputs_stable_during_read") is not True:
        raise ValueError("Input stability must be verified before plotting")
    source_groups = summary["group_counts"]
    group_map = {}
    for group in source_groups:
        displayed = display_group(group)
        if displayed in group_map:
            raise ValueError("Ambiguous group aliases cannot be pooled")
        group_map[displayed] = group
    if set(group_map) != {"Diabetes", "Prediabetes", "Control"}:
        raise ValueError("Exactly Diabetes, Prediabetes and Control groups are required")
    leads = read_table(run_dir / "lead_loci.tsv")
    if not leads or len({r["locus_id"] for r in leads}) != len(leads):
        raise ValueError("A nonempty, unique selected-locus table is required")
    if len(leads) != summary["selected_loci"] or dict(Counter(r["status"] for r in leads)) != summary["locus_status_counts"]:
        raise ValueError("Locus statuses do not match the completed-run summary")
    by_lead = {r["locus_id"]: r for r in leads}
    pairs, frequencies = {}, {}
    haplotypes = {r["locus_id"]: set() for r in leads}
    for source in read_table(run_dir / "haplotype_comparisons.tsv"):
        ident = source["locus_id"]
        if ident not in by_lead or source["group_a"] not in source_groups or source["group_b"] not in source_groups:
            raise ValueError("Pairwise table refers to an unknown locus or group")
        if source["group_a"] == source["group_b"]:
            raise ValueError("A pairwise comparison needs two different groups")
        if by_lead[ident]["status"] != "analysed":
            raise ValueError("An excluded locus cannot contain haplotype comparisons")
        key = ident, source["haplotype_id"], frozenset([source["group_a"], source["group_b"]])
        if key in pairs:
            raise ValueError("Duplicate pairwise comparison, including reversed duplicates")
        pairs[key] = source
        haplotypes[ident].add(source["haplotype_id"])
    for source in read_table(run_dir / "haplotype_group_frequencies.tsv"):
        ident = source["locus_id"]
        if ident not in by_lead or source["group"] not in source_groups:
            raise ValueError("Frequency table refers to an unknown locus or group")
        if by_lead[ident]["status"] != "analysed":
            raise ValueError("An excluded locus cannot contain haplotype frequencies")
        key = ident, source["haplotype_id"], source["group"]
        if key in frequencies:
            raise ValueError("Duplicate haplotype frequency")
        frequencies[key] = numeric(source["haplotype_frequency"], "haplotype_frequency", missing=True)
        if frequencies[key] is not None and not 0 <= frequencies[key] <= 1:
            raise ValueError("Haplotype frequencies must be within [0,1]")
        haplotypes[ident].add(source["haplotype_id"])

    output = []
    for lead in leads:
        ident = lead["locus_id"]
        for hap in sorted(haplotypes[ident]) or [None]:
            for label in ["Diabetes", "Prediabetes"]:
                case, control = group_map[label], group_map["Control"]
                row = dict.fromkeys(FIELDS)
                row.update(locus_id=ident, chrom=lead["chrom"], lead_pos=lead["lead_pos"], haplotype_id=hap,
                           contrast=f"{label} vs Control", case_group=case, control_group=control,
                           panel_status=lead["status"], comparison_status="not_estimable")
                row["haplotype_frequency_case_descriptive"] = frequencies.get((ident, hap, case))
                row["haplotype_frequency_control_descriptive"] = frequencies.get((ident, hap, control))
                source = pairs.get((ident, hap, frozenset([case, control])))
                if lead["status"] != "analysed":
                    row["reason"] = lead["status"]
                elif hap is None:
                    row["reason"] = "no_resolved_haplotypes"
                elif source is None:
                    row["reason"] = "pairwise_result_unavailable"
                else:
                    row["reason"] = source["test_status"]
                    reverse = source["group_a"] == control
                    for destination, origin in [("case", "b" if reverse else "a"), ("control", "a" if reverse else "b")]:
                        for prefix in ["n_resolved", "carriers", "carrier_frequency"]:
                            row[f"{prefix}_{destination}"] = numeric(source[f"{prefix}_{origin}"], f"{prefix}_{origin}", missing=True)
                    if source["test_status"] == "tested_unadjusted":
                        params = summary["parameters"]
                        if not params["test_associations"] or not params["reference_genotypes_confirmed"]:
                            raise ValueError("Tested effects require verified reference-call provenance")
                        for suffix, group in [("case", case), ("control", control)]:
                            n, k, frequency = (row[f"{prefix}_{suffix}"] for prefix in ["n_resolved", "carriers", "carrier_frequency"])
                            if n is None or k is None or frequency is None or n != int(n) or k != int(k) or not 0 <= k <= n <= source_groups[group]:
                                raise ValueError("Invalid participant counts in tested contrast")
                            if n < params["min_test_samples"] or n / source_groups[group] < params["min_resolved_fraction"]:
                                raise ValueError("Tested contrast fails participant-resolution requirements")
                            if not math.isclose(frequency, k / n, abs_tol=1e-7):
                                raise ValueError("Carrier frequency disagrees with participant counts")
                        total = row["n_resolved_case"] + row["n_resolved_control"]
                        carriers = row["carriers_case"] + row["carriers_control"]
                        if min(carriers, total - carriers) < params["min_carriers"]:
                            raise ValueError("Tested contrast fails carrier-count requirements")
                        effect = numeric(source["carrier_frequency_difference_a_minus_b"], "carrier difference")
                        low = numeric(source["carrier_difference_CI95_low"], "difference lower interval")
                        high = numeric(source["carrier_difference_CI95_high"], "difference upper interval")
                        if reverse:
                            effect, low, high = -effect, -high, -low
                        if not -1 - 1e-8 <= low <= effect <= high <= 1 + 1e-8:
                            raise ValueError("Invalid carrier-frequency difference interval")
                        if not math.isclose(effect, row["carrier_frequency_case"] - row["carrier_frequency_control"], abs_tol=1e-7):
                            raise ValueError("Carrier difference disagrees with its group frequencies")
                        p = numeric(source["p_fisher_two_sided"], "p-value")
                        q = numeric(source["q_BH_all_pairwise_tests"], "BH q-value")
                        if not 0 <= p <= q + 1e-8 or not 0 <= q <= 1:
                            raise ValueError("Invalid pairwise p/q values")
                        row.update(comparison_status="tested", carrier_frequency_difference=effect, CI95_low=low,
                                   CI95_high=high, p_fisher_two_sided=p, q_BH_all_pairwise_tests=q)
                output.append(row)
    return summary, output


def short_reason(reason, parameters):
    labels = {
        "lead_below_min_call_rate": f"Lead genotypes below {parameters['min_call_rate']:.0%} completeness",
        "palindromic_lead_orientation_unconfirmed": "Palindromic allele orientation unconfirmed",
        "negative_liftover_allele_orientation_unconfirmed": "Allele orientation after liftover unconfirmed",
        "not_tested_reference_provenance": "Reference-call provenance unconfirmed",
        "too_few_resolved_samples": "Too few resolved participants",
        "low_resolved_fraction": "Too few participants have resolved haplotypes",
        "too_few_carriers_or_noncarriers": "Too few carriers or noncarriers",
        "no_carrier_variation": "No carrier variation among resolved participants",
        "not_requested": "Association test not requested",
        "no_resolved_haplotypes": "No resolved haplotypes",
        "pairwise_result_unavailable": "Pairwise result unavailable",
        "fewer_than_two_eligible_markers": "Fewer than two eligible markers",
        "no_unique_PASS_sequence_allele_match_at_lead": "No eligible exact lead-allele match",
        "multiple_matching_lead_records": "Ambiguous lead-allele match",
    }
    return labels.get(reason, reason.replace("_", " "))


def interpretation(summary, rows):
    tested = [r for r in rows if r["comparison_status"] == "tested"]
    significant = [r for r in tested if r["q_BH_all_pairwise_tests"] < .05]
    lines = ["Haplotype contrasts against Control", "",
             "Requested contrasts: Diabetes vs Control and Prediabetes vs Control.",
             f"{len(tested)} requested haplotype comparisons were tested; {len(rows) - len(tested)} are not estimable."]
    if not tested:
        lines += ["No haplotype difference can be established for either requested contrast from this run.",
                  "No candidate causal variant can be localized from these unavailable haplotype comparisons.",
                  "Not estimable does not mean zero effect, identical groups, or absence of causal variants."]
    else:
        lines += [f"{len(significant)} requested comparisons have BH q < 0.05 in the original all-pairwise family."]
        for r in significant:
            direction = "higher" if r["carrier_frequency_difference"] > 0 else "lower"
            lines.append(f"{r['locus_id']} {r['haplotype_id']}: {r['contrast']}, {direction} carrier frequency; "
                         f"difference {100*r['carrier_frequency_difference']:.2f} percentage points "
                         f"(95% CI {100*r['CI95_low']:.2f} to {100*r['CI95_high']:.2f}), q={r['q_BH_all_pairwise_tests']:.4g}.")
        lines += ["These exploratory associations do not establish causal variants or disease progression."]
    if any(r["comparison_status"] != "tested" for r in rows):
        reasons = Counter(r["reason"] for r in rows if r["comparison_status"] != "tested")
        lines += ["", "Unavailable comparison reasons:"] + [f"{n}: {reason}" for reason, n in sorted(reasons.items())]
    lines += ["", "The effect is participant carrier frequency in Diabetes/Prediabetes minus Control, not chromosome-copy haplotype frequency.",
              "Intervals use pointwise Newcombe-Wilson 95% confidence limits. BH q-values retain the full original pairwise family; they are not recalculated for this picture.",
              "Chromosome-copy haplotype frequencies are separate descriptive fields in contrast_data.tsv.",
              "Recover validated reference/alternate genotypes and usable phase before attempting unavailable comparisons. Candidate localization additionally requires independently validated fine-mapping evidence.",
              "Results, figures and this interpretation are private artifacts."]
    return "\n".join(lines) + "\n"


def draw_contrasts(summary, rows):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    keys = list(dict.fromkeys((r["locus_id"], r["haplotype_id"]) for r in rows))
    lookup = {(r["locus_id"], r["haplotype_id"], r["contrast"]): r for r in rows}
    any_tested = any(r["comparison_status"] == "tested" for r in rows)
    figure, axes = plt.subplots(1, 2, figsize=(15, max(6.5, 3.4 + .43 * len(keys))), sharey=True)
    figure.subplots_adjust(left=.24, right=.91 if any_tested else .97, top=.79, bottom=.18,
                           wspace=.30 if any_tested else .10)
    labels = []
    for key in keys:
        row = lookup[(*key, "Diabetes vs Control")]
        text = f"{row['locus_id']}  {row['chrom']}:{int(row['lead_pos']):,}"
        if row["haplotype_id"]:
            text += f"\n{row['haplotype_id']}"
        labels.append(text)
    for ax, contrast, color in zip(axes, ["Diabetes vs Control", "Prediabetes vs Control"], ["#b64b18", "#008666"]):
        ax.set_title(contrast, fontsize=14, weight="bold", color=color, pad=15)
        ax.set_ylim(len(keys)-.5, -.5)
        for i, key in enumerate(keys):
            row = lookup[(*key, contrast)]
            ax.axhspan(i-.47, i+.47, color="#f4f6f8" if i % 2 == 0 else "white", zorder=0)
            if row["comparison_status"] == "tested":
                effect, low, high = [100*row[k] for k in ["carrier_frequency_difference", "CI95_low", "CI95_high"]]
                ax.errorbar(effect, i, xerr=[[effect-low], [high-effect]], fmt="o", color=color,
                            markersize=5, capsize=3, linewidth=1.3, zorder=3)
                ax.text(1.01, i, f"q={row['q_BH_all_pairwise_tests']:.3g}", transform=ax.get_yaxis_transform(),
                        va="center", ha="left", fontsize=8, clip_on=False)
            else:
                ax.text(.5, i, "Not estimable\n" + short_reason(row["reason"], summary["parameters"]),
                        transform=ax.get_yaxis_transform(), ha="center", va="center", fontsize=8.2, color="#775620",
                        bbox={"boxstyle": "round,pad=.25", "facecolor": "#fff5e3", "edgecolor": "none"}, zorder=4)
        if any_tested:
            ax.set_xlim(-105, 105)
            ax.axvline(0, color="#6d7b89", linestyle="--", linewidth=.8, zorder=1)
            ax.set_xticks([-100, -50, 0, 50, 100])
            ax.set_xlabel("Carrier-frequency difference\n(percentage points; group minus Control)", fontsize=9)
            ax.grid(axis="x", color="#e2e7ec", linewidth=.6)
        else:
            ax.set_xlim(0, 1)
            ax.set_xticks([])
        ax.set_yticks(range(len(keys)), labels)
        ax.tick_params(axis="y", labelsize=8.5, length=0, pad=9)
        for spine in ax.spines.values():
            spine.set_visible(False)
    tested = sum(r["comparison_status"] == "tested" for r in rows)
    figure.text(.055, .957, "Haplotype contrasts against Control", fontsize=21, weight="bold", color="#17334b")
    groups = summary["group_counts"]
    counts = " · ".join(f"{display_group(g)} n={n}" for g, n in sorted(groups.items(), key=lambda item: ["Diabetes", "Prediabetes", "Control"].index(display_group(item[0]))))
    figure.text(.055, .920, counts, fontsize=11.5, color="#526579")
    figure.text(.055, .879, f"{tested} tested comparisons · {len(rows)-tested} not estimable · {summary['selected_loci']} selected regions",
                fontsize=11, color="#526579")
    note = ("Points show participant-carrier differences with pointwise 95% intervals; q-values retain the original pairwise BH family."
            if any_tested else "No effect estimates or significance marks are plotted. Each cell records why that comparison is unavailable.")
    figure.text(.055, .090, note, fontsize=9.5, color="#40566d")
    figure.text(.055, .061, "Unavailable comparisons do not imply identical groups, zero effects, or absence of causal variants.", fontsize=9.5, color="#526579")
    scope = "SYNTHETIC SOFTWARE FIXTURE" if summary.get("synthetic") else "Cohort analysis"
    figure.text(.055, .027, f"PRIVATE RESULTS · {scope} · Carrier tests use participants; chromosome-copy haplotype frequencies are descriptive.",
                fontsize=8.5, color="#63758a")
    return figure


def export_contrasts(run_dir, out):
    run_dir, out = Path(run_dir), Path(out)
    summary, rows = load_contrasts(run_dir)
    if out.exists() and any(out.iterdir()):
        raise ValueError("Use a new or empty contrast output directory")
    source_hashes = {name: hashlib.sha256((run_dir / name).read_bytes()).hexdigest() for name in SOURCES}
    out.mkdir(parents=True, exist_ok=True)
    with (out / "contrast_data.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({key: "NA" if value is None else value for key, value in row.items()} for row in rows)
    figure = draw_contrasts(summary, rows)
    for suffix in ["png", "pdf", "svg"]:
        figure.savefig(out / f"haplotype_contrasts.{suffix}", dpi=300, facecolor="white")
    from matplotlib import pyplot as plt
    plt.close(figure)
    (out / "INTERPRETATION.txt").write_text(interpretation(summary, rows))
    if source_hashes != {name: hashlib.sha256((run_dir / name).read_bytes()).hexdigest() for name in SOURCES}:
        raise ValueError("Source analysis changed during contrast rendering")
    manifest = {"private_results": True, "source_run": run_dir.name, "source_sha256": source_hashes,
                "comparisons_tested": sum(r["comparison_status"] == "tested" for r in rows),
                "comparisons_not_estimable": sum(r["comparison_status"] != "tested" for r in rows),
                "source_completed_at_utc": summary["completed_at_utc"],
                "files": ["haplotype_contrasts.png", "haplotype_contrasts.pdf", "haplotype_contrasts.svg", "contrast_data.tsv", "INTERPRETATION.txt"]}
    (out / "contrast_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True, help="New private local results directory")
    args = parser.parse_args()
    print(json.dumps(export_contrasts(args.run_dir, args.out), indent=2))
