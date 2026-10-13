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
          "panel_status", "comparison_status", "inference_status", "reason", "n_resolved_case", "n_resolved_control",
          "carriers_case", "carriers_control", "carrier_frequency_case", "carrier_frequency_control",
          "descriptive_carrier_frequency_difference",
          "carrier_frequency_difference", "CI95_low", "CI95_high", "p_fisher_two_sided",
          "q_BH_all_pairwise_tests", "haplotype_frequency_case_descriptive", "haplotype_frequency_control_descriptive"]
SOURCES = ["run_summary.json", "lead_loci.tsv", "haplotype_comparisons.tsv", "haplotype_group_frequencies.tsv", "COMPLETED.txt"]
INFERENCE_SOURCE_FIELDS = ["carrier_frequency_difference_a_minus_b", "carrier_difference_CI95_low",
                           "carrier_difference_CI95_high", "p_fisher_two_sided", "q_BH_all_pairwise_tests",
                           "conditional_odds_ratio_a_vs_b", "CI95_low", "CI95_high", "q_BH_all_haplotype_tests"]


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
    if any(isinstance(n, bool) or not isinstance(n, int) or n <= 0 for n in source_groups.values()):
        raise ValueError("Group sizes must be positive integer participant counts")
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
        if source.get("haplotype_id") in {None, "", "NA", "."}:
            raise ValueError("Missing haplotype ID")
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
        if source.get("haplotype_id") in {None, "", "NA", "."}:
            raise ValueError("Missing haplotype ID")
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
                           panel_status=lead["status"], comparison_status="not_estimable", inference_status="not_tested")
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
                    row["inference_status"] = source["test_status"]
                    reverse = source["group_a"] == control
                    for destination, origin in [("case", "b" if reverse else "a"), ("control", "a" if reverse else "b")]:
                        for prefix in ["n_resolved", "carriers", "carrier_frequency"]:
                            row[f"{prefix}_{destination}"] = numeric(source[f"{prefix}_{origin}"], f"{prefix}_{origin}", missing=True)
                    for suffix, group in [("case", case), ("control", control)]:
                        n, k, frequency = (row[f"{prefix}_{suffix}"] for prefix in ["n_resolved", "carriers", "carrier_frequency"])
                        if n is None and k is None and frequency is None:
                            continue
                        if n is None or k is None or n != int(n) or k != int(k) or not 0 <= k <= n <= source_groups[group]:
                            raise ValueError("Invalid participant counts in contrast")
                        if n == 0:
                            if frequency is not None:
                                raise ValueError("Zero resolved participants require an unavailable carrier frequency")
                        elif frequency is None or not math.isclose(frequency, k / n, abs_tol=1e-7):
                            raise ValueError("Carrier frequency disagrees with participant counts")
                        else:
                            # Calculate the descriptive quantity from counts, not rounded source frequencies.
                            row[f"carrier_frequency_{suffix}"] = k / n
                    if all(row[f"carrier_frequency_{suffix}"] is not None for suffix in ["case", "control"]):
                        row["descriptive_carrier_frequency_difference"] = row["carrier_frequency_case"] - row["carrier_frequency_control"]
                        row["comparison_status"] = "descriptive_only"
                    if source["test_status"] != "tested_unadjusted" and any(source.get(field) not in {None, "", "NA", "."} for field in INFERENCE_SOURCE_FIELDS):
                        raise ValueError("Untested comparisons must not contain inferential estimates, intervals or p/q values")
                    if source["test_status"] == "tested_unadjusted":
                        params = summary["parameters"]
                        if not params["test_associations"] or not params["reference_genotypes_confirmed"]:
                            raise ValueError("Tested effects require verified reference-call provenance")
                        for suffix, group in [("case", case), ("control", control)]:
                            n, k, frequency = (row[f"{prefix}_{suffix}"] for prefix in ["n_resolved", "carriers", "carrier_frequency"])
                            if n is None or k is None or frequency is None or n <= 0:
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
    observed = [r for r in rows if r["descriptive_carrier_frequency_difference"] is not None]
    descriptive = [r for r in rows if r["comparison_status"] == "descriptive_only"]
    missing = [r for r in rows if r["comparison_status"] == "not_estimable"]
    significant = [r for r in tested if r["q_BH_all_pairwise_tests"] < .05]
    regions = len({r["locus_id"] for r in observed})
    lines = ["Haplotype contrasts against Control", "",
             "Requested contrasts: Diabetes vs Control and Prediabetes vs Control.",
             f"{len(tested)} requested haplotype comparisons were tested; {len(descriptive)} have descriptive observations only; {len(missing)} lack a computable carrier contrast.",
             f"Observed carrier differences are available in {regions} regions. The complete source analysis selected {summary['selected_loci']} regions.",
             "An observed difference describes the resolved participants only; missing or reference-depleted calls can bias both its numerator and denominator."]
    if not tested:
        lines += ["No statistically supported haplotype association can be established for either requested contrast from this run.",
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
    for label in ["Diabetes vs Control", "Prediabetes vs Control"]:
        candidates = [r for r in observed if r["contrast"] == label]
        if candidates:
            r = max(candidates, key=lambda item: abs(item["descriptive_carrier_frequency_difference"]))
            lines.append(f"Largest absolute observed difference for {label}: {r['haplotype_id']}, "
                         f"{r['carriers_case']:.0f}/{r['n_resolved_case']:.0f} versus {r['carriers_control']:.0f}/{r['n_resolved_control']:.0f} carriers, "
                         f"{100*r['descriptive_carrier_frequency_difference']:+.2f} percentage points. This magnitude ranking is descriptive, not evidence of significance or causality.")
    if descriptive or missing:
        reasons = Counter(r["reason"] for r in descriptive + missing)
        lines += ["", "Reasons inference is unavailable (first failing gate; other limitations may also apply):"]
        lines += [f"{n}: {reason}" for reason, n in sorted(reasons.items())]
    lines += ["", "Observed carrier differences for every available pattern (group minus Control):"]
    for r in observed:
        delta = 100*r["descriptive_carrier_frequency_difference"]
        direction = "higher" if delta > 0 else "lower" if delta < 0 else "the same observed frequency"
        lines.append(f"{r['haplotype_id']} | {r['contrast']} | "
                     f"{r['carriers_case']:.0f}/{r['n_resolved_case']:.0f} versus {r['carriers_control']:.0f}/{r['n_resolved_control']:.0f} | "
                     f"{delta:+.2f} percentage points ({direction}) | {r['comparison_status']} | {r['reason']}")
    lines += ["", "The displayed quantity is participant carrier frequency in Diabetes/Prediabetes minus Control, not chromosome-copy haplotype frequency.",
              "Descriptive-only points have no confidence intervals, p-values, q-values or significance marks. Actual observed zero differences are retained as zero; unavailable comparisons remain missing.",
              "When valid tests exist, intervals use pointwise Newcombe-Wilson 95% confidence limits. BH q-values retain the full original pairwise family; they are not recalculated for this picture.",
              "Chromosome-copy haplotype frequencies are separate descriptive fields in contrast_data.tsv.",
              "Recover validated reference/alternate genotypes and usable phase before attempting unavailable comparisons. Never substitute reference alleles for missing genotypes. Candidate localization additionally requires independently validated fine-mapping evidence.",
              "Results, figures and this interpretation are private artifacts."]
    return "\n".join(lines) + "\n"


def _contrast_value(row):
    """Return a plotted value only for an explicitly available observation."""
    if row["comparison_status"] == "tested":
        return row.get("carrier_frequency_difference")
    if row["comparison_status"] == "descriptive_only":
        return row.get("descriptive_carrier_frequency_difference")
    return None


def _draw_contrast_rows(summary, rows, *, regional=False):
    """Two comparison columns with counts and honest inference availability."""
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    keys = list(dict.fromkeys((r["locus_id"], r["haplotype_id"]) for r in rows))
    if not keys:
        raise ValueError("No contrast rows are available to draw")
    lookup = {(r["locus_id"], r["haplotype_id"], r["contrast"]): r for r in rows}
    observed = sum(_contrast_value(r) is not None for r in rows)
    tested = sum(r["comparison_status"] == "tested" for r in rows)
    height = max(7.0, 3.8 + .50 * len(keys))
    figure, axes = plt.subplots(1, 2, figsize=(18 if observed else 15, height), sharey=True)
    figure.subplots_adjust(left=.14 if regional else .23, right=.985, top=.77, bottom=.22, wspace=.12)
    labels = []
    for key in keys:
        row = lookup[(*key, "Diabetes vs Control")]
        if regional and row["haplotype_id"]:
            label = row["haplotype_id"].removeprefix(row["locus_id"] + "_")
        else:
            label = f"{row['locus_id']}  {row['chrom']}:{int(row['lead_pos']):,}"
            if row["haplotype_id"]:
                label += f"\n{row['haplotype_id']}"
        labels.append(label)
    for ax, contrast, color in zip(axes, ["Diabetes vs Control", "Prediabetes vs Control"], ["#b64b18", "#008666"]):
        ax.set_title(contrast, fontsize=13, weight="bold", color=color, pad=23 if observed else 15)
        ax.set_ylim(len(keys)-.5, -.5)
        if observed:
            ax.text(111, 1.02, "Carriers / resolved participants\nGroup vs Control; difference", transform=ax.get_xaxis_transform(),
                    ha="left", va="bottom", fontsize=8, color="#526579")
        for i, key in enumerate(keys):
            row = lookup[(*key, contrast)]
            ax.axhspan(i-.47, i+.47, color="#f4f6f8" if i % 2 == 0 else "white", zorder=0)
            value = _contrast_value(row)
            if row["comparison_status"] == "tested":
                effect, low, high = [100*row[k] for k in ["carrier_frequency_difference", "CI95_low", "CI95_high"]]
                ax.errorbar(effect, i, xerr=[[effect-low], [high-effect]], fmt="o", color=color,
                            markersize=5, capsize=3, linewidth=1.3, zorder=3)
            elif value is not None:
                ax.plot(100*value, i, marker="D", color=color, markersize=5, linestyle="none", zorder=3)
            if value is not None:
                counts = f"{int(row['carriers_case'])}/{int(row['n_resolved_case'])} vs {int(row['carriers_control'])}/{int(row['n_resolved_control'])}"
                label = f"{counts}\n{100*value:+.2f} pp; "
                label += f"q={row['q_BH_all_pairwise_tests']:.3g}" if row["comparison_status"] == "tested" else "observed only"
                ax.text(111, i, label, ha="left", va="center", fontsize=8.1, color="#233445")
            else:
                ax.text(.5, i, "Not estimable\n" + short_reason(row["reason"], summary["parameters"]),
                        transform=ax.get_yaxis_transform(), ha="center", va="center", fontsize=8.2, color="#775620",
                        bbox={"boxstyle": "round,pad=.25", "facecolor": "#fff5e3", "edgecolor": "none"}, zorder=4)
        if observed:
            ax.set_xlim(-105, 205)
            ax.axvline(0, color="#6d7b89", linestyle="--", linewidth=.8, zorder=1)
            ax.set_xticks([-100, -50, 0, 50, 100])
            ax.set_xlabel("Carrier-frequency difference (percentage points; group minus Control)", fontsize=8.5)
            for tick in [-100, -50, 50, 100]:
                ax.axvline(tick, color="#e2e7ec", linewidth=.6, zorder=0)
        else:
            ax.set_xlim(0, 1)
            ax.set_xticks([])
        ax.set_yticks(range(len(keys)), labels)
        ax.tick_params(axis="y", labelsize=8.5, length=0, pad=9)
        for spine in ax.spines.values():
            spine.set_visible(False)
    title = "Haplotype contrasts against Control"
    if regional:
        first = rows[0]
        title = f"{first['locus_id']} · {first['chrom']}:{int(first['lead_pos']):,} · Haplotype carrier contrasts"
    figure.text(.045, .957, title, fontsize=19, weight="bold", color="#17334b")
    groups = summary["group_counts"]
    counts = " · ".join(f"{display_group(g)} n={n}" for g, n in sorted(groups.items(), key=lambda item: ["Diabetes", "Prediabetes", "Control"].index(display_group(item[0]))))
    figure.text(.045, .913, counts, fontsize=11, color="#526579")
    scope = "DESCRIPTIVE ONLY" if observed and not tested else "INFERENCE UNAVAILABLE" if not tested else "OBSERVED AND TESTED CONTRASTS"
    figure.text(.045, .866, f"{scope} · {observed} observed contrasts · {tested} tested · {len(rows)-observed} unavailable", fontsize=11, color="#526579")
    if observed:
        note = "Diamonds: observed carrier differences among resolved participants; no confidence intervals or p-values are implied."
        if tested:
            note += " Circles: tested differences with 95% intervals."
    else:
        note = "No observed haplotype differences are available for these rows. Each cell records the reason."
    figure.text(.045, .123, note, fontsize=9, color="#40566d")
    if observed and not tested:
        reasons = sorted({short_reason(r["reason"], summary["parameters"]) for r in rows if _contrast_value(r) is not None})
        reason_note = "; ".join(reasons)
        if len(reason_note) > 155:
            reason_note = "Reasons are recorded per comparison in contrast_data.tsv."
        figure.text(.045, .088, "Inference unavailable: " + reason_note, fontsize=9, color="#775620")
    else:
        figure.text(.045, .088, "Unavailable comparisons do not imply identical groups, zero effects, or absence of causal variants.", fontsize=9, color="#526579")
    scope = "SYNTHETIC SOFTWARE FIXTURE" if summary.get("synthetic") else "Cohort analysis"
    figure.text(.045, .051, "Observed values condition on resolved calls and may be biased by selective missingness. No causal conclusion follows.", fontsize=9, color="#526579")
    figure.text(.045, .023, f"PRIVATE RESULTS · {scope} · Carrier counts use participants, not their two chromosome copies.", fontsize=8.5, color="#63758a")
    return figure


def draw_region_contrasts(summary, rows):
    """Return one regional forest figure with counts, descriptive points and valid inference."""
    if len({r["locus_id"] for r in rows}) != 1:
        raise ValueError("A regional contrast picture requires exactly one region")
    return _draw_contrast_rows(summary, rows, regional=True)


def _draw_unavailable_summary(summary, rows):
    """Summarize large unavailable runs without constructing an enormous bitmap."""
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    import textwrap
    reasons = sorted({row["reason"] for row in rows})
    counts = Counter((row["reason"], row["contrast"]) for row in rows)
    figure, axes = plt.subplots(1, 2, figsize=(15, max(6.5, 3.5 + .55*len(reasons))), sharey=True)
    figure.subplots_adjust(left=.36, right=.97, bottom=.22, top=.73, wspace=.15)
    for ax, contrast in zip(axes, ["Diabetes vs Control", "Prediabetes vs Control"]):
        ax.set_title(contrast, fontsize=13, weight="bold", color="#17334b")
        ax.set_xlim(0, 1); ax.set_ylim(len(reasons)-.5, -.5); ax.set_xticks([])
        labels = [textwrap.fill(short_reason(reason, summary["parameters"]), width=45) for reason in reasons]
        ax.set_yticks(range(len(reasons)), labels, fontsize=9)
        ax.tick_params(length=0)
        for i, reason in enumerate(reasons):
            ax.axhspan(i-.46, i+.46, color="#f4f6f8", zorder=0)
            ax.text(.5, i, f"{counts[reason, contrast]} unavailable rows", ha="center", va="center", fontsize=10, color="#775620")
        for spine in ax.spines.values():
            spine.set_visible(False)
    nregions = len({row["locus_id"] for row in rows})
    figure.text(.045, .95, "Haplotype contrasts: observed values unavailable", fontsize=20, weight="bold", color="#17334b")
    figure.text(.045, .875, f"{nregions} regions · {len(rows)} unavailable rows · Exclusion counts shown below", fontsize=11, color="#526579")
    figure.text(.045, .12, "No observed differences or inferential estimates are drawn. Every region and its reason remain in contrast_data.tsv.", fontsize=9.5, color="#526579")
    figure.text(.045, .075, "These are counts of unavailable rows, not effects. Unavailable does not mean zero difference or equal groups.", fontsize=9.5, color="#526579")
    scope = "SYNTHETIC SOFTWARE FIXTURE" if summary.get("synthetic") else "Cohort analysis"
    figure.text(.045, .03, "PRIVATE RESULTS · " + scope, fontsize=9, color="#63758a")
    return figure


def draw_contrasts(summary, rows):
    """Return a scalable overview, separating observed differences from inference."""
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.colors import Normalize
    import numpy as np
    observed = [r for r in rows if _contrast_value(r) is not None]
    keys = list(dict.fromkeys((r["locus_id"], r["haplotype_id"]) for r in rows))
    if not observed and len(keys) > 50:
        return _draw_unavailable_summary(summary, rows)
    if not observed or len(keys) <= 12:
        return _draw_contrast_rows(summary, rows)
    loci = list(dict.fromkeys(r["locus_id"] for r in observed))
    included = [r for r in rows if r["locus_id"] in loci and r["haplotype_id"] is not None]
    lookup = {(r["locus_id"], r["haplotype_id"], r["contrast"]): r for r in included}
    haps = {locus: list(dict.fromkeys(r["haplotype_id"] for r in included if r["locus_id"] == locus)) for locus in loci}
    cols = min(3, len(loci)); nrows = math.ceil(len(loci)/cols)
    heights = [1.05 + .31*max(len(haps[locus]) for locus in loci[index*cols:(index+1)*cols]) for index in range(nrows)]
    height = sum(heights) + 3.8
    figure = plt.figure(figsize=(17.5, height))
    grid = figure.add_gridspec(nrows, cols, left=.065, right=.985, top=1-2.55/height, bottom=1.6/height,
                              hspace=.50, wspace=.37, height_ratios=heights)
    maximum = max(abs(100*_contrast_value(r)) for r in observed)
    bound = min(100, max(5, 5*math.ceil(maximum/5)))
    cmap = plt.get_cmap("RdBu_r").with_extremes(bad="#cfd5dc")
    norm = Normalize(-bound, bound)
    for index, locus in enumerate(loci):
        ax = figure.add_subplot(grid[index//cols, index%cols])
        patterns = haps[locus]
        cells = [[lookup[locus, hap, contrast] for contrast in ["Diabetes vs Control", "Prediabetes vs Control"]] for hap in patterns]
        values = np.asarray([[np.nan if _contrast_value(row) is None else 100*_contrast_value(row) for row in pair] for pair in cells])
        ax.imshow(np.ma.masked_invalid(values), cmap=cmap, norm=norm, interpolation="nearest", aspect="auto")
        for y, pair in enumerate(cells):
            for x, row in enumerate(pair):
                value = _contrast_value(row)
                label = "NA" if value is None else f"{100*value:+.2f}" if value != 0 else "0.00"
                if row["comparison_status"] == "tested":
                    label += " †"
                color = "white" if value is not None and abs(100*value) > .65*bound else "#233445"
                ax.text(x, y, label, ha="center", va="center", fontsize=8.5, color=color)
        ax.set_xticks([0, 1], ["Diabetes − Control\npercentage points", "Prediabetes − Control\npercentage points"], fontsize=8.5)
        ax.set_yticks(range(len(patterns)), [hap.removeprefix(locus + "_") for hap in patterns], fontsize=8.5)
        ax.tick_params(length=0)
        first = cells[0][0]
        ax.set_title(f"{locus} · {first['chrom']}:{int(first['lead_pos']):,}\n{len(patterns)} local haplotypes", fontsize=10.5, weight="bold", pad=10)
        for spine in ax.spines.values():
            spine.set_visible(False)
    tested = sum(r["comparison_status"] == "tested" for r in rows)
    excluded = len({r["locus_id"] for r in rows}) - len(loci)
    label = "DESCRIPTIVE ONLY" if not tested else "OBSERVED CARRIER DIFFERENCES"
    synthetic = "SYNTHETIC FIXTURE · " if summary.get("synthetic") else ""
    figure.text(.045, 1-.32/height, synthetic + "Haplotype carrier differences versus Control", fontsize=20, weight="bold", color="#17334b", va="top")
    figure.text(.045, 1-.77/height, f"{label} · {len(loci)} regions with observations · {sum(len(haps[locus]) for locus in loci)} local haplotypes · {len(observed)} observed contrasts · {tested} tested", fontsize=11, color="#526579", va="top")
    cax = figure.add_axes([.30, 1-1.55/height, .43, .14/height])
    colorbar = figure.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
    colorbar.set_label("Observed participant-carrier frequency difference (percentage points; group − Control)", fontsize=8.5)
    colorbar.ax.tick_params(labelsize=8)
    figure.text(.045, .88/height, f"{excluded} additional regions have no observed contrasts; their exclusion reasons remain in contrast_data.tsv. Gray = unavailable; 0.00 = observed zero.", fontsize=9, color="#526579")
    figure.text(.045, .62/height, "Each value is carriers / resolved participants in the disease group minus Control. Counts and denominators appear in separate regional pictures.", fontsize=9, color="#526579")
    if tested:
        note = "† Tested comparisons: see regional pictures for 95% confidence intervals and original BH q-values. Other cells are descriptive only."
    else:
        note = "No inferential tests, confidence intervals or adjusted p-values are available. Test-blocking reasons appear in contrast_data.tsv."
    figure.text(.045, .36/height, note, fontsize=9, color="#775620")
    figure.text(.045, .10/height, "PRIVATE RESULTS · Selective missingness can bias observed values. Local haplotype IDs are region-specific; no causal conclusion follows.", fontsize=9, color="#63758a")
    return figure


def export_contrasts(run_dir, out):
    run_dir, out = Path(run_dir), Path(out)
    source_hashes = {name: hashlib.sha256((run_dir / name).read_bytes()).hexdigest() for name in SOURCES}
    summary, rows = load_contrasts(run_dir)
    if out.exists() and any(out.iterdir()):
        raise ValueError("Use a new or empty contrast output directory")
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
    index = ["# Haplotype carrier contrasts by region", "",
             "Observed differences are descriptive unless explicitly supported by a valid test. Gray/unavailable values are not zero.", "",
             "| Region | Haplotypes observed | Figure or exclusion |", "| --- | ---: | --- |"]
    regional_files = []
    for ident in dict.fromkeys(row["locus_id"] for row in rows):
        regional = [r for r in rows if r["locus_id"] == ident]
        available = [r for r in regional if r["descriptive_carrier_frequency_difference"] is not None]
        if not available:
            index.append(f"| {ident} | 0 | {short_reason(regional[0]['reason'], summary['parameters'])} |")
            continue
        folder = out / "regional" / ident
        folder.mkdir(parents=True)
        regional_figure = draw_region_contrasts(summary, regional)
        for suffix in ["png", "pdf"]:
            path = folder / f"haplotype_contrasts.{suffix}"
            regional_figure.savefig(path, dpi=220, facecolor="white")
            regional_files.append(str(path.relative_to(out)))
        plt.close(regional_figure)
        (folder / "INTERPRETATION.txt").write_text(interpretation(summary, regional))
        regional_files.append(str((folder / "INTERPRETATION.txt").relative_to(out)))
        n_haplotypes = len({r['haplotype_id'] for r in available})
        index.append(f"| [{ident}](regional/{ident}/INTERPRETATION.txt) | {n_haplotypes} | [PNG](regional/{ident}/haplotype_contrasts.png) · [PDF](regional/{ident}/haplotype_contrasts.pdf) |")
    (out / "index.md").write_text("\n".join(index) + "\n")
    if source_hashes != {name: hashlib.sha256((run_dir / name).read_bytes()).hexdigest() for name in SOURCES}:
        raise ValueError("Source analysis changed during contrast rendering")
    manifest = {"private_results": True, "source_run": run_dir.name, "source_sha256": source_hashes,
                "comparisons_tested": sum(r["comparison_status"] == "tested" for r in rows),
                "comparisons_descriptive_only": sum(r["comparison_status"] == "descriptive_only" for r in rows),
                "comparisons_not_estimable": sum(r["comparison_status"] == "not_estimable" for r in rows),
                "selected_regions": summary["selected_loci"],
                "regions_with_observed_contrasts": len(regional_files) // 3,
                "source_completed_at_utc": summary["completed_at_utc"],
                "files": ["haplotype_contrasts.png", "haplotype_contrasts.pdf", "haplotype_contrasts.svg", "contrast_data.tsv", "INTERPRETATION.txt", "index.md"] + regional_files}
    (out / "contrast_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True, help="New private local results directory")
    args = parser.parse_args()
    print(json.dumps(export_contrasts(args.run_dir, args.out), indent=2))
