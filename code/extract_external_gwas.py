#!/usr/bin/env python3
"""Extract unfiltered SNP/small-indel statistics from a declared external GWAS.

The supported profile is the GRCh38 UKB TOPMed PheWeb T2D download (250.2).
This step preserves supplied REF/ALT alleles and text precision; it does not
normalize variants, calculate LD, fit a model, or change cohort genotypes.
"""
from __future__ import annotations

import argparse
from bisect import bisect_right
from collections import Counter, OrderedDict
import csv
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
import gzip
import hashlib
import json
import math
from pathlib import Path
import platform
import re
import shutil


PROFILE = "ukb-topmed-t2d-250.2"
FIELDS = ["locus_id", "variant_id", "chrom", "pos", "ref", "alt", "variant_type",
          "beta", "se", "pval", "af", "case_af", "control_af", "rsids",
          "nearest_genes", "consequence", "source_line"]
SOURCE_FIELDS = {"chrom", "pos", "ref", "alt", "beta", "sebeta", "pval", "af",
                 "case_af", "control_af", "rsids", "nearest_genes", "consequence"}
WINDOW_FIELDS = ["locus_id", "chrom", "lead_pos", "start", "end"]


def fingerprint(path):
    path = Path(path)
    before = path.stat()
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError("Input changed while hashing")
    return dict(path=str(path.resolve()), bytes=after.st_size,
                mtime_ns=after.st_mtime_ns, sha256=digest.hexdigest())


def open_text(path):
    with Path(path).open("rb") as stream:
        compressed = stream.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt", newline="") if compressed else Path(path).open(newline="")


def classify(ref, alt, limit=49):
    """Classify stored explicit sequences without changing representation."""
    if not ref or not alt or not set(ref + alt) <= set("ACGT") or ref == alt:
        return "non_sequence_or_identical_alleles"
    if len(ref) == len(alt) == 1:
        return "SNP"
    if len(ref) == len(alt):
        return "MNV"
    return "small_indel" if max(len(ref), len(alt)) <= limit + 1 else "larger_or_complex_indel"


def read_windows(path):
    """Require safe, unique, nonoverlapping 1-based inclusive autosomal windows."""
    with Path(path).open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)) or not set(WINDOW_FIELDS) <= set(fields):
            raise ValueError("Lead table needs unique locus_id/chrom/lead_pos/start/end columns")
        rows = list(reader)
    if not rows:
        raise ValueError("A nonempty lead-window table is required")
    identifiers, regions = set(), {}
    for row in rows:
        if None in row or any(row.get(k) is None for k in fields):
            raise ValueError("Malformed lead-window row")
        ident = row["locus_id"]
        if not re.fullmatch(r"[A-Za-z0-9_-]+", ident) or ident in identifiers:
            raise ValueError("Safe, unique locus identifiers are required")
        identifiers.add(ident)
        chrom = row["chrom"]
        if chrom not in {f"chr{i}" for i in range(1, 23)}:
            raise ValueError("Lead chromosomes must use chr1 through chr22")
        try:
            start, end, lead = (int(row[k]) for k in ["start", "end", "lead_pos"])
        except (TypeError, ValueError) as exc:
            raise ValueError("Window coordinates must be positive integers") from exc
        if not 1 <= start <= lead <= end:
            raise ValueError("Each positive inclusive window must contain its lead")
        regions.setdefault(chrom, []).append((start, end, ident))
    for windows in regions.values():
        windows.sort()
        if any(left[1] >= right[0] for left, right in zip(windows, windows[1:])):
            raise ValueError("Extraction requires nonoverlapping source windows")
    return rows, regions, {chrom: [r[0] for r in windows] for chrom, windows in regions.items()}


def valid_statistics(row):
    """Do not underflow tiny signed p-values through a float conversion."""
    try:
        beta, se = float(row["beta"]), float(row["sebeta"])
        if not math.isfinite(beta) or not math.isfinite(se) or se <= 0:
            return False
        for name in ["pval", "af", "case_af", "control_af"]:
            number = Decimal(row[name])
            if not number.is_finite() or not 0 <= number <= 1:
                return False
    except (TypeError, ValueError, InvalidOperation):
        return False
    return True


class RegionWriters:
    """Bound file descriptors independently of the number of requested windows."""
    def __init__(self, folder, identifiers, limit=32):
        self.folder, self.limit, self.opened = folder, limit, OrderedDict()
        for ident in identifiers:
            with (folder / f"{ident}.tsv").open("w", newline="") as stream:
                csv.DictWriter(stream, FIELDS, delimiter="\t", lineterminator="\n").writeheader()

    def write(self, ident, row):
        if ident not in self.opened:
            if len(self.opened) >= self.limit:
                self.opened.popitem(last=False)[1][0].close()
            stream = (self.folder / f"{ident}.tsv").open("a", newline="")
            self.opened[ident] = stream, csv.DictWriter(stream, FIELDS, delimiter="\t", lineterminator="\n")
        self.opened.move_to_end(ident)
        self.opened[ident][1].writerow(row)

    def close(self):
        for stream, _writer in self.opened.values():
            stream.close()
        self.opened.clear()


def source_profile(name):
    if name != PROFILE:
        raise ValueError("Unsupported external GWAS source profile")
    cases, controls = 18105, 387087
    return dict(
        profile=name, download_url="https://pheweb.org/UKB-TOPMed/download/250.2",
        source_about_url="https://pheweb.org/UKB-TOPMed/about",
        phenotype_url="https://pheweb.org/UKB-TOPMed/pheno/250.2",
        source_description="UKB TOPMed-imputed White British participants, SAIGE T2D phecode 250.2",
        source_build="GRCh38", effect_allele="ALT", effect_scale="log_odds_ratio",
        case_count_from_portal=cases, control_count_from_portal=controls,
        total_participants_from_portal=cases + controls,
        n_eff_balanced_case_control_formula=4 / (1 / cases + 1 / controls),
        n_eff_note="Trait-level case/control effective size; not measured variant-specific sample size",
        portal_metadata_verified_on="2026-10-13",
        covariates="Relatedness, sex, birth year, first four principal components",
        genotype_availability="External summary statistics only; not cohort genotypes",
        source_identity_user_declared=True, byte_identity_to_remote_release_verified=False,
        provenance_note="The required source-profile option declares the local file's identity. Schema checks and a local SHA-256 do not attest identity to the remote release.")


def extract(gwas, leads, out, *, source=PROFILE, max_indel_length=49, download_provenance=None):
    gwas, leads, out = (Path(p).resolve() for p in [gwas, leads, out])
    if out.exists():
        raise ValueError("Use a new output directory; previous outputs are never overwritten")
    if not isinstance(max_indel_length, int) or not 1 <= max_indel_length <= 49:
        raise ValueError("max_indel_length must be an integer from 1 to 49")
    profile = source_profile(source)
    before = {"gwas": fingerprint(gwas), "leads": fingerprint(leads)}
    loci, regions, starts = read_windows(leads)
    if download_provenance is not None:
        before["download_provenance"] = fingerprint(download_provenance)
    # Validate schema before creating outputs; incomplete scans have no completion marker.
    with open_text(gwas) as stream:
        fields = next(csv.reader(stream, delimiter="\t"), [])
    if len(fields) != len(set(fields)) or not SOURCE_FIELDS <= set(fields):
        raise ValueError("Unexpected UKB-TOPMed source header")
    out.mkdir(parents=True)
    shutil.copyfile(leads, out / "lead_loci.tsv")
    folder = out / "regions"
    folder.mkdir()
    writers = RegionWriters(folder, [row["locus_id"] for row in loci])
    region_counts = {row["locus_id"]: Counter() for row in loci}
    counts, global_types, seen = Counter(), Counter(), set()
    try:
        with open_text(gwas) as stream, gzip.open(out / "regional_variants.tsv.gz", "wt", newline="") as dest:
            writer = csv.DictWriter(dest, FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for line, row in enumerate(csv.DictReader(stream, delimiter="\t"), 2):
                if None in row or any(row.get(k) is None for k in fields):
                    raise ValueError(f"Malformed source row at line {line}")
                counts["source_rows"] += 1
                chrom = row["chrom"].removeprefix("chr")
                try:
                    pos = int(row["pos"])
                except ValueError as exc:
                    raise ValueError(f"Invalid source position at line {line}") from exc
                if pos < 1 or not re.fullmatch(r"(?:[1-9]|1[0-9]|2[0-2]|X|Y|MT|M)", chrom):
                    raise ValueError(f"Invalid source coordinate at line {line}")
                chrom = "chr" + chrom
                kind = classify(row["ref"], row["alt"], max_indel_length)
                global_types[kind] += 1
                idx = bisect_right(starts.get(chrom, []), pos) - 1
                if idx < 0 or pos > regions[chrom][idx][1]:
                    continue
                ident = regions[chrom][idx][2]
                counter = region_counts[ident]
                counts["rows_in_requested_windows"] += 1
                counter["source_rows"] += 1
                counter["source_" + kind] += 1
                if kind not in {"SNP", "small_indel"}:
                    counts["excluded_" + kind] += 1
                    continue
                variant = f"{chrom}:{pos}:{row['ref']}:{row['alt']}"
                if variant in seen:
                    raise ValueError(f"Duplicate exact source variant {variant}; resolve before fine-mapping")
                seen.add(variant)
                if not valid_statistics(row):
                    counts["excluded_invalid_summary_statistics"] += 1
                    counter["excluded_invalid_summary_statistics"] += 1
                    continue
                output = dict(locus_id=ident, variant_id=variant, chrom=chrom, pos=pos,
                    ref=row["ref"], alt=row["alt"], variant_type=kind, beta=row["beta"],
                    se=row["sebeta"], pval=row["pval"], af=row["af"], case_af=row["case_af"],
                    control_af=row["control_af"], rsids=row["rsids"], nearest_genes=row["nearest_genes"],
                    consequence=row["consequence"], source_line=line)
                writers.write(ident, output)
                writer.writerow(output)
                counter["included_" + kind] += 1
                counts["included_" + kind] += 1
    finally:
        writers.close()
    for record in before.values():
        stat = Path(record["path"]).stat()
        if (stat.st_size, stat.st_mtime_ns) != (record["bytes"], record["mtime_ns"]):
            raise ValueError("Input changed during extraction; discard incomplete output")
    if fingerprint(out / "lead_loci.tsv")["sha256"] != before["leads"]["sha256"]:
        raise ValueError("Copied lead table differs from source provenance")
    counter_fields = ["source_rows", "source_SNP", "source_small_indel", "source_MNV",
                      "source_non_sequence_or_identical_alleles", "source_larger_or_complex_indel",
                      "included_SNP", "included_small_indel", "excluded_invalid_summary_statistics"]
    with (out / "region_summary.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, WINDOW_FIELDS + counter_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in loci:
            writer.writerow({**{k: row[k] for k in WINDOW_FIELDS},
                             **{k: region_counts[row["locus_id"]][k] for k in counter_fields}})
    manifest = dict(private_results=True, created_at_utc=datetime.now(timezone.utc).isoformat(),
        analysis="external_summary_statistics_extraction_not_a_finemapping_fit", python_version=platform.python_version(),
        source_files=before, script=fingerprint(__file__), source_provenance=profile,
        source_allele_class_counts=dict(global_types), extraction_counts=dict(counts),
        selected_regions=len(loci), region_files_written=len(loci), max_indel_length=max_indel_length,
        window_coordinate_system="1-based inclusive; nonoverlapping autosomal windows",
        indel_gate="Unequal ACGT allele sequence lengths; each stored allele length <= max_indel_length+1",
        p_value_filter_applied=False, reference_alleles_normalized=False,
        ld_available_in_this_artifact=False, new_statistical_tests_or_fits=0,
        limitations=[
            "Only requested windows are extracted; these need not cover complete fine-mapping loci.",
            "Allele representations and rounded source beta/SE values remain unchanged; normalization and LD/z diagnostics are required.",
            "Source-profile declaration is not byte-level verification of the local file's download provenance.",
            "External GWAS estimates do not establish Diabetes/Prediabetes/Control contrasts in the local cohort.",
            "Trait-level case/control counts are not measured variant-specific sample sizes.",
            "Any subsequent credible sets are conditional on retained variants, LD, models and windows."],
        output_files={name: fingerprint(out / name) for name in
                      ["lead_loci.tsv", "regional_variants.tsv.gz", "region_summary.tsv"]},
        region_files={row["locus_id"]: fingerprint(folder / f"{row['locus_id']}.tsv") for row in loci})
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    (out / "COMPLETED.txt").write_text("External statistics extracted; no new fit performed. Read manifest.json.\n")
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["gwas", "leads", "out"]:
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--source-profile", choices=[PROFILE], required=True,
                        help="Declare and independently verify the local download's source identity")
    parser.add_argument("--max-indel-length", type=int, default=49)
    parser.add_argument("--download-provenance", type=Path,
                        help="Optional local download/provenance record to fingerprint; contents are not copied")
    args = parser.parse_args()
    manifest = extract(args.gwas, args.leads, args.out, source=args.source_profile,
                       max_indel_length=args.max_indel_length, download_provenance=args.download_provenance)
    print(json.dumps({k: manifest[k] for k in ["selected_regions", "extraction_counts"]}, indent=2))


if __name__ == "__main__":
    main()
