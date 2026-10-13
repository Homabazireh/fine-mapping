#!/usr/bin/env python3
"""Audit saved fine-mapping evidence without converting unvalidated PIPs to claims.

Only aggregate counts and public GWAS rows are exported. This is an evidence
audit, not a new fine-mapping fit. No legacy model is eligible for causal ranking
without independent validation of its LD, harmonization and source provenance.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path


METHOD_SOURCES = [
    {"title": "SuSiE-RSS summary-statistic and LD diagnostics",
     "url": "https://stephenslab.github.io/susieR/articles/susierss_diagnostic.html",
     "application": "Validate that GWAS allele coding and reference LD are consistent before interpreting fine-mapping."},
    {"title": "SuSiE-RSS reference documentation",
     "url": "https://stephenslab.github.io/susieR/reference/susie_rss.html",
     "application": "GWAS sample size is recommended, although optional; omission alone does not invalidate a model."},
]

LEGACY_FIELDS = [
    "legacy_locus", "run_status", "qc_present", "model_converged", "ld_sample_size",
    "gwas_n", "ld_shrinkage", "variants_analyzed", "pip_rows_present",
    "credible_set_memberships_present", "retained_genotype_rows_found",
    "retained_gt_ref_hom", "retained_gt_het", "retained_gt_alt_hom",
    "retained_gt_missing_or_other", "retained_gt_width_mismatch_rows",
    "allele_or_z_alignment_errors", "palindromic_retained_variants",
    "eligible_for_causal_ranking", "exclusion_reasons",
]
LEAD_FIELDS = [
    "locus_id", "chrom", "position", "effect_allele", "other_allele",
    "external_gwas_p", "external_gwas_beta", "external_gwas_se", "external_gwas_eaf",
    "external_gwas_ncases", "external_gwas_ncontrols", "external_gwas_neff",
    "effect_direction", "liftover_strand", "allele_orientation_review",
    "cohort_locus_status", "cohort_association_evidence", "causal_status",
]
CANDIDATE_FIELDS = ["locus_id", "variant", "evidence", "validated_pip", "causal_status"]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def fingerprint(path: Path) -> dict:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return {"name": path.name, "bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def normalized_chrom(value: str) -> str:
    return "chr" + value.removeprefix("chr")


def external_leads(leads: list[dict], gwas_file: Path) -> list[dict]:
    """Join exact ordered EA/NEA pairs and retain tiny p-values as strings."""
    wanted = {(normalized_chrom(r["chrom"]), r["lead_pos"],
               r["lead_effect_allele"].upper(), r["lead_other_allele"].upper()) for r in leads}
    found = {}
    with gwas_file.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Chromsome", "Position", "EffectAllele", "NonEffectAllele", "Pval", "Beta", "SE"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("GWAS audit requires the project's original GWAS column names")
        for row in reader:
            key = (normalized_chrom(row["Chromsome"]), row["Position"],
                   row["EffectAllele"].upper(), row["NonEffectAllele"].upper())
            if key in wanted:
                if key in found:
                    raise ValueError("Ambiguous duplicate GWAS record for a selected lead")
                found[key] = row
    result = []
    for lead in leads:
        key = (normalized_chrom(lead["chrom"]), lead["lead_pos"],
               lead["lead_effect_allele"].upper(), lead["lead_other_allele"].upper())
        if key not in found:
            raise ValueError(f"Selected lead lacks exact GWAS allele match: {lead['locus_id']}")
        source = found[key]
        # String comparison is deliberately avoided: Decimal retains ultra-small P.
        if Decimal(source["Pval"]) != Decimal(lead["lead_p"]):
            raise ValueError("Lead/GWAS p-value disagreement: source files may differ")
        beta, se = float(source["Beta"]), float(source["SE"])
        if not math.isfinite(beta) or not math.isfinite(se) or se <= 0:
            raise ValueError("Non-finite GWAS effect or non-positive SE")
        palindromic = set(key[2:]) in ({"A", "T"}, {"C", "G"})
        strand = source.get("liftover_strand", "NA")
        orientation = ("palindromic_requires_forward_reference_validation" if palindromic
                       else "negative_liftover_strand_requires_harmonization" if strand == "-"
                       else "exact_allele_match_only_not_reference_attestation")
        result.append({
            "locus_id": lead["locus_id"], "chrom": key[0], "position": key[1],
            "effect_allele": key[2], "other_allele": key[3],
            "external_gwas_p": source["Pval"], "external_gwas_beta": source["Beta"],
            "external_gwas_se": source["SE"], "external_gwas_eaf": source.get("EAF", "NA"),
            "external_gwas_ncases": source.get("Ncases", "NA"),
            "external_gwas_ncontrols": source.get("Ncontrols", "NA"),
            "external_gwas_neff": source.get("Neff", "NA"),
            "effect_direction": "positive" if beta > 0 else "negative" if beta < 0 else "zero",
            "liftover_strand": strand, "allele_orientation_review": orientation,
            "cohort_locus_status": lead["status"],
            "cohort_association_evidence": "not_estimated",
            "causal_status": "unresolved_external_association_only",
        })
    return result


def alignment_errors(rows: list[dict]) -> tuple[int, int]:
    errors = palindromic = 0
    for row in rows:
        ea, nea = row.get("GWAS_EA", "").upper(), row.get("GWAS_NEA", "").upper()
        ref, alt = row.get("ONT_REF", "").upper(), row.get("ONT_ALT", "").upper()
        palindromic += set((ref, alt)) in ({"A", "T"}, {"C", "G"})
        try:
            beta, se, z = (float(row[name]) for name in ("Beta", "SE", "Z_ALT"))
            expected = beta / se * (1 if ea == alt else -1)
            ok = (set((ea, nea)) == set((ref, alt)) and ref != alt
                  and se > 0 and math.isfinite(expected) and math.isfinite(z)
                  and math.isclose(z, expected, rel_tol=1e-6, abs_tol=1e-8))
        except (ValueError, KeyError, ZeroDivisionError):
            ok = False
        errors += not ok
    return errors, palindromic


def retained_genotype_counts(path: Path, keys: set[str], expected_n: int | None) -> dict:
    counts = Counter()
    if not path.exists():
        return {"available": False}
    seen = set()
    with path.open() as handle:
        for line in handle:
            parts = line.rstrip("\n\r").split("\t")
            if len(parts) < 5:
                raise ValueError("Malformed saved legacy genotype row")
            chrom, pos, ref, alt = parts[:4]
            a1, a2 = sorted((ref.upper(), alt.upper()))
            key = f"{chrom}:{pos}:{a1}:{a2}"
            if key not in keys:
                continue
            if key in seen:
                raise ValueError("Duplicate retained legacy genotype variant")
            seen.add(key)
            counts["rows"] += 1
            counts["width_mismatch_rows"] += expected_n is not None and len(parts) - 4 != expected_n
            for gt, n in Counter(parts[4:]).items():
                canonical = gt.replace("|", "/")
                kind = {"0/0": "ref_hom", "0/1": "het", "1/0": "het", "1/1": "alt_hom"}.get(
                    canonical, "missing_or_other")
                counts[kind] += n
    return {"available": True, **{k: counts[k] for k in
            ["rows", "ref_hom", "het", "alt_hom", "missing_or_other", "width_mismatch_rows"]}}


def audit_legacy(legacy_dir: Path) -> tuple[list[dict], dict]:
    summary_path = legacy_dir / "all_loci_summary.tsv"
    if not summary_path.is_file():
        return [], {"available": False, "reason": "legacy_summary_not_supplied"}
    script = legacy_dir / "run_one_locus.R"
    script_text = script.read_text() if script.exists() else ""
    method = {
        "available": True, "summary_source": fingerprint(summary_path),
        "script_source": fingerprint(script) if script.exists() else None,
        "script_snapshot_is_not_proof_of_historical_execution": True,
        "script_mean_imputes_missing_dosages": "dosage[i,missing] <- mean(" in script_text,
        "script_has_dp_gq_filters": "min_dp" in script_text or "min_gq" in script_text,
        "script_explicitly_checks_ld_z_consistency": any(
            token in script_text for token in ("estimate_s_rss(", "kriging_rss(")),
        "script_mentions_liftover_strand_or_palindromic_checks": any(
            token in script_text for token in ("liftover_strand", "palindrom")),
    }
    result = []
    for record in read_tsv(summary_path):
        name = f"{record['CHR']}_{record['LEAD_BP']}"
        locus = legacy_dir / "loci" / name
        qc_rows = read_tsv(locus / "susie_QC.tsv") if (locus / "susie_QC.tsv").exists() else []
        if len(qc_rows) > 1:
            raise ValueError("Expected one legacy QC row per locus")
        qc = qc_rows[0] if qc_rows else {}
        pip_rows = read_tsv(locus / "susie_PIP.tsv") if (locus / "susie_PIP.tsv").exists() else []
        keys = {r["Variant_key"] for r in pip_rows}
        if len(keys) != len(pip_rows):
            raise ValueError("Duplicate legacy PIP variant keys")
        errors, palindromic = alignment_errors(pip_rows)
        n = int(qc["ONT_samples"]) if qc.get("ONT_samples", "").isdigit() else None
        counts = retained_genotype_counts(locus / "locus_genotypes.tsv", keys, n)
        cs = read_tsv(locus / "susie_credible_sets.tsv") if (locus / "susie_credible_sets.tsv").exists() else []
        reasons = ["ld_callability_and_gwas_population_match_unvalidated",
                   "ld_z_consistency_not_validated", "functional_causal_evidence_not_supplied"]
        if record["STATUS"] != "SUCCESS":
            reasons.append("legacy_fit_failed")
        if not qc:
            reasons.append("legacy_qc_missing")
        elif qc.get("Converged") != "TRUE":
            reasons.append("model_not_converged")
        if qc.get("GWAS_n") == "NOT_SUPPLIED":
            reasons.append("gwas_sample_size_not_supplied")
        if keys and not counts["available"]:
            reasons.append("ld_genotype_artifact_missing")
        elif keys:
            if counts["rows"] != len(keys):
                reasons.append("retained_ld_variants_missing_from_genotypes")
            if counts["width_mismatch_rows"]:
                reasons.append("ld_sample_count_mismatch")
            if counts["ref_hom"] == 0:
                reasons.append("no_reference_homozygotes_in_retained_ld_genotypes")
        if errors:
            reasons.append("allele_or_z_alignment_error")
        if palindromic:
            reasons.append("palindromic_orientation_unvalidated")
        result.append({
            "legacy_locus": name, "run_status": record["STATUS"], "qc_present": bool(qc),
            "model_converged": qc.get("Converged", "NA"), "ld_sample_size": n if n is not None else "NA",
            "gwas_n": qc.get("GWAS_n", "NA"), "ld_shrinkage": qc.get("LD_shrinkage", "NA"),
            "variants_analyzed": qc.get("Variants_analyzed", "NA"), "pip_rows_present": len(pip_rows),
            "credible_set_memberships_present": len(cs),
            "retained_genotype_rows_found": counts.get("rows", "NA"),
            **{f"retained_gt_{k}": counts.get(k, "NA") for k in
               ["ref_hom", "het", "alt_hom", "missing_or_other", "width_mismatch_rows"]},
            "allele_or_z_alignment_errors": errors, "palindromic_retained_variants": palindromic,
            "eligible_for_causal_ranking": False, "exclusion_reasons": ";".join(reasons),
        })
    return result, method


def audit_finemapping(run_dir: Path, gwas: Path, legacy_dir: Path, out: Path) -> dict:
    """Write an aggregate audit for the current run with no tested haplotypes."""
    run_dir, gwas, legacy_dir, out = map(Path, (run_dir, gwas, legacy_dir, out))
    if out.exists() and any(out.iterdir()):
        raise ValueError("Use a new or empty output folder")
    summary_path = run_dir / "run_summary.json"
    source = json.loads(summary_path.read_text())
    if not (run_dir / "COMPLETED.txt").exists():
        raise ValueError("Cohort run is incomplete")
    if source.get("synthetic") or "demo" in run_dir.name.lower():
        raise ValueError("Synthetic results cannot supply cohort causal evidence")
    leads = read_tsv(run_dir / "lead_loci.tsv")
    if source.get("haplotype_carrier_tests_performed", 0) or any(
            int(r.get("selected_markers", "0")) > 0 for r in leads):
        raise ValueError("This audit is for unresolved cohorts; evaluate tested haplotypes separately")
    evidence = external_leads(leads, gwas)
    legacy, method = audit_legacy(legacy_dir)
    old_haplotype_script = legacy_dir.parent / "Haplotype_comparison" / "run_one_locus_haplotype.R"
    historical_haplotype = {"available": old_haplotype_script.is_file()}
    if old_haplotype_script.is_file():
        historical_text = old_haplotype_script.read_text()
        historical_haplotype.update({
            "source": fingerprint(old_haplotype_script),
            "effect_allele_assigned_risk_without_beta_sign": (
                'risk_code_anchor <- if (g$Effect_allele_in_ONT[anchor_i] == "ALT") 1 else 0'
                in historical_text and "g$Beta" not in historical_text),
            "gwas_p_converted_to_machine_numeric": "as.numeric(g$Pval)" in historical_text,
            "phase_set_continuity_validation_present": "phase_set" in historical_text,
            "results_used_in_current_analysis": False,
        })
    reasons = Counter(reason for row in legacy for reason in row["exclusion_reasons"].split(";"))
    result = {
        "schema_version": 1, "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "analysis": "aggregate_evidence_audit_not_new_finemapping", "synthetic": False,
        "genome_build_declared_in_cohort_run": source.get("genome_build_user_declared", "not_supplied"),
        "cohort_group_counts": source.get("group_counts", {}),
        "cohort_selected_loci": len(leads), "cohort_haplotype_association_tests": 0,
        "cohort_locus_status_counts": dict(Counter(r["status"] for r in leads)),
        "cohort_high_quality_reference_homozygotes": source.get("scan_counts", {}).get(
            "high_quality_homozygous_reference_genotypes"),
        "external_gwas_leads_reported": len(evidence), "resolved_candidate_variants": 0,
        "candidate_status": "not_resolvable_from_validated_available_evidence",
        "legacy_loci": len(legacy), "legacy_status_counts": dict(Counter(r["run_status"] for r in legacy)),
        "legacy_eligible_for_causal_ranking": 0, "legacy_exclusion_counts": dict(reasons),
        "legacy_method": method, "legacy_haplotype_script_review": historical_haplotype,
        "method_sources": METHOD_SOURCES,
        "provenance": {"cohort_summary": fingerprint(summary_path), "gwas": fingerprint(gwas),
                       "selected_leads": fingerprint(run_dir / "lead_loci.tsv")},
        "limitations": [
            "No eligible haplotypes or cohort association tests: no haplotype can be declared different between groups.",
            "External GWAS leads are associations, not evidence of a diabetes/prediabetes/control difference in this cohort.",
            "Legacy convergence or high PIP cannot resolve unvalidated LD and allele provenance; numeric legacy PIPs are withheld.",
            "Absent reference homozygotes warrants callability investigation; it does not prove every missing genotype is reference.",
            ("Missing genotypes remain missing in this audit; the saved legacy script mean-imputes dosage for LD, not missing-to-reference filling."
             if method.get("script_mean_imputes_missing_dosages") else
             "Missing genotypes remain missing in this audit; historical missing-call handling needs source validation."),
            "A missing GWAS sample size is a limitation, not by itself proof of invalid fine-mapping.",
            "The reported effect direction is the sign of the supplied Beta; phenotype coding and effect scale require source verification.",
            "GRCh38 coordinate provenance does not establish genome-wide forward-reference allele orientation.",
            "Legacy variant and genotype counts are locus memberships; overlapping loci may count the same variant more than once.",
        ],
        "required_for_candidate_resolution": [
            "Recover validated reference and alternate calls and continuous phase for cohort haplotypes; retain missingness QC.",
            "Use complete regional summary statistics, verify source trait/effect scale, per-variant sample sizes and GRCh38 harmonization.",
            "Use ancestry-matched validated LD, check allele order and palindromic/strand ambiguity, and evaluate LD/z consistency.",
            "Rerun fine-mapping with recorded settings, sensitivity analyses and calibrated credible sets after those validations.",
            "Integrate functional or colocalization evidence before making stronger causal claims; statistical candidates remain candidates.",
        ],
    }
    if historical_haplotype.get("effect_allele_assigned_risk_without_beta_sign"):
        result["limitations"].append(
            "Historical haplotype code labels the GWAS effect allele as risk without using Beta sign, "
            "converts tiny p-values to machine numeric, and does not validate phase-set continuity. "
            "Its group comparisons are not reused as current haplotype evidence.")
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / "external_gwas_leads.tsv", evidence, LEAD_FIELDS)
    write_tsv(out / "legacy_locus_audit.tsv", legacy, LEGACY_FIELDS)
    write_tsv(out / "candidate_variants.tsv", [], CANDIDATE_FIELDS)
    result["output_sha256"] = {name: fingerprint(out / name)["sha256"] for name in
                               ["external_gwas_leads.tsv", "legacy_locus_audit.tsv", "candidate_variants.tsv"]}
    (out / "audit_summary.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    report = ["Variant-evidence audit", "", "No causal variants can be resolved from the validated available evidence.",
              f"Cohort: {len(leads)} selected loci; 0 eligible haplotype panels; 0 haplotype association tests.",
              f"External GWAS: {len(evidence)} selected lead associations, reported separately from cohort evidence.",
              f"Legacy SuSiE: {len(legacy)} loci audited; 0 eligible for causal ranking.",
              "candidate_variants.tsv is intentionally header-only. This means unresolved, not absence of causal variants.",
              "", "Limitations:", *["- " + value for value in result["limitations"]],
              "", "Evidence needed:", *["- " + value for value in result["required_for_candidate_resolution"]],
              "", "Method sources:", *[f"- {s['title']}: {s['url']}" for s in METHOD_SOURCES]]
    (out / "REPORT.txt").write_text("\n".join(report) + "\n")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--gwas", type=Path, required=True)
    parser.add_argument("--legacy-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    summary = audit_finemapping(args.run_dir, args.gwas, args.legacy_dir, args.out)
    print(f"Audited {summary['legacy_loci']} legacy loci; causal candidate resolution remains unavailable.")


if __name__ == "__main__":
    main()
