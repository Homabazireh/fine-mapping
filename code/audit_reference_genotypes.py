#!/usr/bin/env python3
"""Independently audit raw reference/missing calls without participant output.

Uses indexed bcftools queries and streams raw GTs to awk for bounded-memory
counting. No association gate, genotype, or input file is changed. Output is
aggregate research data and belongs only in a private/local results directory.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile


GT_KEYS = (
    "all_GT_cells", "missing_GT", "fully_missing_GT", "partially_missing_GT",
    "raw_called_GT", "ref_homozygotes", "alt_homozygotes", "heterozygotes",
    "non_diploid_GT", "phased_called", "unphased_called",
)


def fingerprint(path: Path) -> dict:
    before = path.stat()
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise RuntimeError("An input changed while being fingerprinted")
    return {"name": path.name, "bytes": after.st_size,
            "mtime_ns": after.st_mtime_ns, "sha256": digest.hexdigest()}


def classify_gt_counts(genotypes: dict[str, int]) -> dict[str, int]:
    """A partial genotype is missing, never a reference or complete call."""
    counts = Counter({key: 0 for key in GT_KEYS})
    for genotype, n in genotypes.items():
        if n < 0 or not re.fullmatch(r"(?:\d+|\.)(?:[/|](?:\d+|\.))*", genotype):
            raise ValueError("Invalid raw genotype or count")
        alleles = re.split(r"[/|]", genotype)
        counts["all_GT_cells"] += n
        if "." in alleles:
            counts["missing_GT"] += n
            counts["fully_missing_GT" if set(alleles) == {"."}
                   else "partially_missing_GT"] += n
            continue
        counts["raw_called_GT"] += n
        if len(alleles) != 2:
            counts["non_diploid_GT"] += n
        elif set(alleles) == {"0"}:
            counts["ref_homozygotes"] += n
        elif len(set(alleles)) == 1:
            counts["alt_homozygotes"] += n
        else:
            counts["heterozygotes"] += n
        counts["phased_called" if "|" in genotype else "unphased_called"] += n
    return dict(counts)


def merge_windows(leads: list[dict]) -> list[tuple[str, int, int]]:
    regions = defaultdict(list)
    for row in leads:
        start, end = int(row["start"]), int(row["end"])
        if start < 1 or end < start:
            raise ValueError("Region coordinates must be 1-based inclusive")
        regions[row["chrom"]].append((start, end))
    output = []
    for chrom in sorted(regions):
        merged = []
        for start, end in sorted(regions[chrom]):
            if merged and start <= merged[-1][1] + 1:
                merged[-1][1] = max(end, merged[-1][1])
            else:
                merged.append([start, end])
        output.extend((chrom, start, end) for start, end in merged)
    return output


def sanitized_header(header: str) -> tuple[list[str], int]:
    """Allowlist provenance; never save #CHROM or ##SAMPLE identifiers."""
    prefixes = ("##source=", "##reference=", "##longphaseVersion=",
                "##bcftools_mergeVersion=", "##bcftools_mergeCommand=")
    provenance = [line for line in header.splitlines() if line.startswith(prefixes)]
    columns = [line.split("\t") for line in header.splitlines()
               if line.startswith("#CHROM\t")]
    if len(columns) != 1 or len(columns[0]) < 10:
        raise ValueError("A multisample VCF header is required")
    return provenance, len(columns[0]) - 9


def summarize_lead_lines(lines, leads: list[dict], min_dp: float, min_gq: float,
                         n_samples: int) -> dict:
    wanted = defaultdict(set)
    for row in leads:
        wanted[(row["chrom"], row["lead_pos"])].add(
            frozenset((row["lead_effect_allele"], row["lead_other_allele"])))
    genotypes = Counter()
    counts = Counter({"matched_records": 0, "called_DP_missing_or_lt_threshold": 0,
                      "called_DP_pass_GQ_missing_or_lt_threshold": 0,
                      "called_DP_GQ_pass": 0})
    missing = Counter()
    completeness = []
    for line in lines:
        fields = line.rstrip("\n").split("\t")
        chrom, pos, ref, alt, filt = fields[:5]
        # Exact PASS biallelic allele-pair matching, without strand inference.
        if filt != "PASS" or "," in alt or frozenset((ref, alt)) not in wanted[(chrom, pos)]:
            continue
        if len(fields) - 5 != n_samples:
            raise ValueError("Raw lead query sample width differs from header")
        counts["matched_records"] += 1
        called = high_quality = 0
        for cell in fields[5:]:
            genotype, dp, gq, _phase_set = cell.split(":")
            genotypes[genotype] += 1
            if "." in re.split(r"[/|]", genotype):
                missing[f"missing_GT_DP_{'absent' if dp == '.' else 'present'}"
                        f"_GQ_{'absent' if gq == '.' else 'present'}"] += 1
                continue
            called += 1
            if dp == "." or float(dp) < min_dp:
                counts["called_DP_missing_or_lt_threshold"] += 1
            elif gq == "." or float(gq) < min_gq:
                counts["called_DP_pass_GQ_missing_or_lt_threshold"] += 1
            else:
                counts["called_DP_GQ_pass"] += 1
                high_quality += 1
        completeness.append((called, high_quality))
    return {"counts": dict(counts), "raw_GT_counts": classify_gt_counts(genotypes),
            "missing_field_counts": dict(missing),
            "record_completeness": {
                "min_raw_called": min((n for n, _ in completeness), default=None),
                "max_raw_called": max((n for n, _ in completeness), default=None),
                "records_all_samples_raw_called": sum(n == n_samples for n, _ in completeness),
                "records_at_least_90_percent_raw_called": sum(n / n_samples >= .9 for n, _ in completeness),
                "records_at_least_90_percent_DP_GQ_pass": sum(n / n_samples >= .9 for _, n in completeness),
            }}


def write_regions(path: Path, regions) -> None:
    path.write_text("".join(f"{chrom}\t{start}\t{end}\n" for chrom, start, end in regions))


def run_audit(vcf: Path, leads_path: Path, out: Path, min_dp=10., min_gq=20.) -> dict:
    if out.exists():
        raise ValueError("Use a new output directory to preserve earlier audits")
    for program in ("bcftools", "awk"):
        if not shutil.which(program):
            raise RuntimeError(f"Required program is unavailable: {program}")
    with leads_path.open(newline="") as handle:
        leads = list(csv.DictReader(handle, delimiter="\t"))
    if not leads:
        raise ValueError("The lead table is empty")
    inputs = {"vcf": fingerprint(vcf), "lead_loci": fingerprint(leads_path)}
    provenance, n_samples = sanitized_header(subprocess.check_output(
        ["bcftools", "view", "--no-version", "-h", str(vcf)], text=True))
    version = subprocess.check_output(["bcftools", "--version"], text=True).splitlines()[:2]
    windows = merge_windows(leads)
    with tempfile.TemporaryDirectory(prefix="reference-genotype-audit-") as directory:
        window_path = Path(directory) / "windows.tsv"
        lead_path = Path(directory) / "lead_positions.tsv"
        write_regions(window_path, windows)
        write_regions(lead_path, sorted({(r["chrom"], int(r["lead_pos"]), int(r["lead_pos"])) for r in leads}))
        query = subprocess.Popen(["bcftools", "query", "--regions-overlap", "0",
            "-R", str(window_path), "-f", "[%GT\n]", str(vcf)], stdout=subprocess.PIPE)
        try:
            aggregate = subprocess.run(["awk", "{gt[$0]++} END {for (g in gt) print g,gt[g]}"],
                stdin=query.stdout, stdout=subprocess.PIPE, text=True, check=True)
        finally:
            query.stdout.close()
            if query.wait() != 0:
                raise RuntimeError("The raw genotype query failed")
        raw = {g: int(n) for g, n in (line.split() for line in aggregate.stdout.splitlines())}
        all_counts = classify_gt_counts(raw)
        if all_counts["all_GT_cells"] % n_samples:
            raise RuntimeError("Raw genotype cell count does not match VCF sample width")
        query = subprocess.Popen(["bcftools", "query", "--regions-overlap", "0", "-R", str(lead_path),
            "-f", "%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%GT:%DP:%GQ:%PS]\n", str(vcf)],
            stdout=subprocess.PIPE, text=True)
        try:
            lead_counts = summarize_lead_lines(query.stdout, leads, min_dp, min_gq, n_samples)
        finally:
            query.stdout.close()
            if query.wait() != 0:
                raise RuntimeError("The raw lead genotype query failed")
    for key, path in (("vcf", vcf), ("lead_loci", leads_path)):
        stat = path.stat()
        if (stat.st_size, stat.st_mtime_ns) != (inputs[key]["bytes"], inputs[key]["mtime_ns"]):
            raise RuntimeError("An input changed during the audit")
    result = {
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "method": "Independent indexed bcftools query of raw GT, before any genotype quality filtering",
        "sample_count": n_samples, "selected_regions": len(leads),
        "union_windows": len(windows), "records_in_union_windows": all_counts["all_GT_cells"] // n_samples,
        "raw_window_GT_counts": all_counts, "exact_PASS_biallelic_leads": lead_counts,
        "lead_DP_GQ_thresholds_for_diagnosis_only": {"min_dp": min_dp, "min_gq": min_gq},
        "inputs": inputs, "script": fingerprint(Path(__file__)),
        "source_header_provenance": provenance, "bcftools_version": version,
        "query_templates": [
            "bcftools query --regions-overlap 0 -R union_windows.tsv -f '[%GT\\n]' INPUT.vcf.gz",
            "awk '{gt[$0]++} END {for (g in gt) print g,gt[g]}'",
            "bcftools query --regions-overlap 0 -R lead_positions.tsv -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER[\\t%GT:%DP:%GQ:%PS]\\n' INPUT.vcf.gz",
        ],
        "interpretation": {
            "reference_homozygotes_observed_before_QC": all_counts["ref_homozygotes"] > 0,
            "missing_calls_reclassified_as_reference": False,
            "original_variant_only_export_proven": False,
            "reference_provenance_validated_by_this_audit": False,
            "limits": "The source header documents a merge, not whether source VCFs contained every confident reference call. Raw GT counts diagnose missing reference information but cannot reconstruct it, establish associations, or identify causal variants.",
        },
    }
    out.mkdir(parents=True)
    (out / "genotype_audit.json").write_text(json.dumps(result, indent=2) + "\n")
    write_regions(out / "union_windows.tsv", windows)
    shutil.copyfile(Path(__file__), out / "audit_reference_genotypes.py")
    counts = all_counts
    report = ["Raw genotype audit", "", f"Regions: {len(leads)}; records: {result['records_in_union_windows']}; samples: {n_samples}.",
              f"Raw genotype cells: {counts['all_GT_cells']}; missing allele(s): {counts['missing_GT']}; fully called: {counts['raw_called_GT']}.",
              f"Reference homozygotes before DP/GQ filtering: {counts['ref_homozygotes']}.", ""]
    if counts["ref_homozygotes"] == 0:
        report.extend(["Reference homozygotes are already absent from the selected VCF windows before quality filtering; genotype quality thresholds did not remove them.",
                       "The absent calls cannot safely be classified as non-carriers. Obtain reference-inclusive gVCFs or regenerate confident genotypes from original alignments before interpreting population haplotype contrasts.", ""])
    report.extend(["Missing genotypes were not changed. No tests or causal rankings were calculated by this audit.",
                   "Source merge provenance supports investigation of reference calls omitted by a variant-only merge, but the original export command is unavailable; that specific cause is not proven.",
                   "All source file SHA-256 digests were freshly calculated. No participant identifiers are exported.", "",
                   "Reproduce from the repository root, using a new private/local output directory:",
                   f"python code/audit_reference_genotypes.py --vcf {vcf} --leads {leads_path} --out PRIVATE_NEW_OUTPUT",
                   ""])
    (out / "REPORT.txt").write_text("\n".join(report))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--leads", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--min-dp", type=float, default=10.)
    parser.add_argument("--min-gq", type=float, default=20.)
    args = parser.parse_args()
    result = run_audit(args.vcf, args.leads, args.out, args.min_dp, args.min_gq)
    print(json.dumps({"out": str(args.out), "regions": result["selected_regions"],
                      "reference_homozygotes": result["raw_window_GT_counts"]["ref_homozygotes"]}))


if __name__ == "__main__":
    main()
