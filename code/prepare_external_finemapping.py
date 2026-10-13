#!/usr/bin/env python3
"""Prepare exact-allele external GWAS/reference LD inputs without fitting models.

Only complete, biallelic reference genotypes contribute to signed ALT-dosage
correlations. Missing calls are never imputed. This external reference workflow
does not alter cohort association eligibility or cohort genotype filters.
"""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess

import numpy as np


MATCH_FIELDS = ["variant_id", "chrom", "pos", "ref", "alt", "variant_type",
                "beta", "se", "pval", "af", "gwas_n", "reference_af", "source_variant_id"]
QC_FIELDS = ["source_variant_id", "normalized_variant_id", "variant_type", "status", "reasons"]


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path, rows, fields):
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def variant_kind(ref, alt):
    if not ref or not alt or not set(ref + alt) <= set("ACGT") or ref == alt:
        return "unsupported"
    if len(ref) == len(alt) == 1:
        return "SNP"
    return "INDEL" if len(ref) != len(alt) else "unsupported"


def variant_key(chrom, pos, ref, alt):
    return chrom, int(pos), ref, alt


def variant_id(key):
    return ":".join(map(str, key))


def common_af(value, minimum=.01):
    try:
        af = float(value)
    except (TypeError, ValueError):
        return False
    return math.isfinite(af) and 0 <= af <= 1 and min(af, 1-af) >= minimum


def complete_dosages(genotypes, expected_samples):
    """Return observed ALT counts only; reject an entire incomplete record."""
    if len(genotypes) != expected_samples:
        raise ValueError("Reference genotype width differs from verified sample order")
    result = []
    for gt in genotypes:
        alleles = gt.replace("|", "/").split("/")
        if len(alleles) != 2 or any(a not in {"0", "1"} for a in alleles):
            return None
        result.append(int(alleles[0]) + int(alleles[1]))
    return np.asarray(result, dtype=float)


def signed_ld(dosages):
    """Variables are rows and reference participants are columns."""
    x = np.asarray(dosages, dtype=float)
    if x.ndim != 2 or min(x.shape) < 2 or not np.isfinite(x).all():
        raise ValueError("At least two finite variants and participants are required")
    if np.any(np.var(x, axis=1) <= 0):
        raise ValueError("Monomorphic reference variants cannot define LD")
    matrix = np.corrcoef(x)
    if not np.isfinite(matrix).all() or not np.allclose(matrix, matrix.T, atol=1e-12):
        raise ValueError("Invalid signed LD matrix")
    if not np.allclose(np.diag(matrix), 1, atol=1e-12) or np.max(np.abs(matrix)) > 1 + 1e-12:
        raise ValueError("Invalid LD diagonal or range")
    # A correlation Gram matrix is PSD by construction; the fitting implementation
    # independently checks eigenvalues to distinguish roundoff from invalid input.
    return matrix


def normalize_gwas(rows, folder, chrom, fasta, contig_length):
    """Use original row IDs; bcftools -c x excludes REF mismatches without flips."""
    source_vcf, normalized = folder / "gwas.source.vcf", folder / "gwas.normalized.vcf"
    with source_vcf.open("w") as stream:
        stream.write(f"##fileformat=VCFv4.2\n##contig=<ID={chrom},length={contig_length}>\n")
        stream.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for number, row in enumerate(rows):
            if variant_kind(row["ref"], row["alt"]) == "unsupported":
                continue
            if row["chrom"] != chrom or int(row["pos"]) < 1:
                raise ValueError("GWAS row lies on an unexpected chromosome or position")
            stream.write(f"{chrom}\t{row['pos']}\tGWAS{number}\t{row['ref']}\t{row['alt']}\t.\tPASS\t.\n")
    command = ["bcftools", "norm", "--no-version", "-f", str(fasta), "-c", "x",
               "-Ov", "-o", str(normalized), str(source_vcf)]
    with (folder / "normalization.log").open("w") as log:
        subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
    normalized_rows = {}
    with normalized.open() as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            number = int(fields[2].removeprefix("GWAS"))
            if number in normalized_rows:
                raise ValueError("Normalization produced multiple records for one GWAS variant")
            normalized_rows[number] = variant_key(fields[0], fields[1], fields[3], fields[4])
    return normalized_rows, command


def read_reference(path, chrom, start, end, sample_count):
    command = ["bcftools", "query", "--regions-overlap", "0", "-r", f"{chrom}:{start}-{end}",
               "-f", "%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%GT]\n", str(path)]
    text = subprocess.check_output(command, text=True)
    records, key_counts = {}, Counter()
    for line in text.splitlines():
        fields = line.split("\t")
        key = variant_key(*fields[:4])
        key_counts[key] += 1
        reasons = []
        if fields[4] not in {"PASS", "."}:
            reasons.append("reference_site_filter")
        if variant_kind(key[2], key[3]) == "unsupported":
            reasons.append("reference_not_biallelic_sequence_variant")
        dosage = complete_dosages(fields[5:], sample_count)
        if dosage is None:
            reasons.append("reference_incomplete_or_nonbiallelic_GT")
        records[key] = (dosage, reasons)
    for key, count in key_counts.items():
        if count > 1:
            records[key] = (None, ["reference_duplicate_normalized_key"])
    return records, command, key_counts


def match_variants(rows, normalized, reference, start, end, nominal_n,
                   min_maf=.01, max_af_difference=.1):
    """Return aligned rows/dosages and every observable blocker per GWAS row."""
    duplicates = Counter(normalized.values())
    matched, qc, selected_dosages = [], [], []
    for number, row in enumerate(rows):
        key = normalized.get(number)
        kind = variant_kind(row["ref"], row["alt"])
        reasons = []
        if kind == "unsupported":
            reasons.append("unsupported_GWAS_allele_class")
        elif key is None:
            reasons.append("GWAS_reference_allele_mismatch")
        if key is not None:
            if not start <= key[1] <= end:
                reasons.append("normalized_position_outside_window")
            if duplicates[key] > 1:
                reasons.append("GWAS_duplicate_normalized_key")
        try:
            beta, se, p, af = (float(row[name]) for name in ["beta", "se", "pval", "af"])
            valid = all(math.isfinite(v) for v in [beta, se, p, af]) and se > 0 and 0 <= p <= 1 and 0 <= af <= 1
        except (ValueError, KeyError):
            valid = False
        if not valid:
            reasons.append("invalid_GWAS_summary_statistics")
        elif min(af, 1-af) < min_maf:
            reasons.append("GWAS_MAF_below_threshold")
        dosage, ref_af = None, None
        if key is not None:
            if key not in reference:
                reasons.append("no_exact_normalized_reference_match")
            else:
                dosage, reference_reasons = reference[key]
                reasons.extend(reference_reasons)
                if dosage is not None:
                    ref_af = float(np.mean(dosage)/2)
                    if float(np.var(dosage)) <= 0:
                        reasons.append("reference_no_dosage_variation")
                    if min(ref_af, 1-ref_af) < min_maf:
                        reasons.append("reference_MAF_below_threshold")
                    if valid and abs(af-ref_af) > max_af_difference + 1e-12:
                        reasons.append("reference_GWAS_AF_difference_above_threshold")
        qc.append(dict(source_variant_id=row["variant_id"], normalized_variant_id=variant_id(key) if key else "NA",
                       variant_type=kind, status="excluded" if reasons else "matched", reasons=";".join(reasons) or "none"))
        if not reasons:
            if dosage is None:
                raise ValueError("A matched variant must have observed reference dosages")
            matched.append(dict(variant_id=variant_id(key), chrom=key[0], pos=key[1], ref=key[2], alt=key[3],
                                variant_type=variant_kind(key[2], key[3]), beta=row["beta"], se=row["se"],
                                pval=row["pval"], af=row["af"], gwas_n=nominal_n, reference_af=ref_af,
                                source_variant_id=row["variant_id"]))
            selected_dosages.append(dosage)
    order = sorted(range(len(matched)), key=lambda i: (matched[i]["pos"], matched[i]["ref"], matched[i]["alt"]))
    return [matched[i] for i in order], [selected_dosages[i] for i in order], qc


def prepare(source_dir, reference_dir, fasta, out, workers=4, locus=None):
    source_dir, reference_dir, fasta, out = map(lambda p: Path(p).resolve(), [source_dir, reference_dir, fasta, out])
    if out.exists():
        raise ValueError("Use a new output directory")
    if not 1 <= workers <= 8:
        raise ValueError("Use 1..8 workers")
    for directory in [source_dir, reference_dir]:
        if not (directory / "COMPLETED.txt").is_file():
            raise ValueError("Only completed source/reference acquisitions may be prepared")
    source = json.loads((source_dir / "manifest.json").read_text())
    ref = json.loads((reference_dir / "manifest.json").read_text())
    input_hashes = {str(path): sha(path) for path in [source_dir/"manifest.json", reference_dir/"manifest.json"]}
    gwas = source["source_provenance"]
    if gwas["source_build"] != "GRCh38" or gwas["effect_allele"] != "ALT" or ref["genome_build"] != "GRCh38":
        raise ValueError("Both sources must declare GRCh38 and GWAS effects must be ALT-oriented")
    if ref.get("population") != "EUR_original_2504_unrelated_subset" or not ref.get("upstream_filtered_release_ACCEPTS_PASS_or_dot"):
        raise ValueError("Unexpected reference population or unconfirmed filtered-release convention")
    if sha(fasta) != ref["fasta_sha256"]:
        raise ValueError("FASTA differs from reference normalization provenance")
    input_hashes[str(fasta)] = ref["fasta_sha256"]
    if sha(source_dir / "lead_loci.tsv") != source["source_files"]["leads"]["sha256"]:
        raise ValueError("Source lead table checksum mismatch")
    panel = reference_dir / "source_sample_panel.tsv"
    if sha(panel) != ref["sample_panel_sha256"]:
        raise ValueError("Reference population metadata checksum mismatch")
    samples = [row["sample"] for row in read_tsv(panel) if row["super_pop"] == "EUR"]
    if len(samples) != 503 or len(set(samples)) != 503 or ref["sample_count"] != 503:
        raise ValueError("Exactly 503 unique original unrelated EUR samples are required")
    if (reference_dir / "EUR_503.samples.txt").read_text().splitlines() != samples:
        raise ValueError("Reference sample list/order differs from population metadata")
    leads = read_tsv(source_dir / "lead_loci.tsv")
    if not leads or len({r["locus_id"] for r in leads}) != len(leads) or any(
            re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", r["locus_id"]) is None for r in leads):
        raise ValueError("Unique nonempty safe locus identifiers are required")
    input_hashes[str(source_dir/"lead_loci.tsv")] = source["source_files"]["leads"]["sha256"]
    input_hashes[str(panel)] = ref["sample_panel_sha256"]
    selected = [row for row in leads if locus is None or row["locus_id"] == locus]
    if not selected:
        raise ValueError("Requested locus is absent")
    with Path(str(fasta)+".fai").open() as stream:
        contig_lengths = {v[0]: int(v[1]) for v in (line.split("\t") for line in stream)}
    references = {entry["chrom"]: entry for entry in ref["chromosomes"]}
    for chrom in {row["chrom"] for row in selected}:
        entry = references[chrom]
        path = reference_dir / entry["normalized_file"]
        if sha(path) != entry["normalized_sha256"]:
            raise ValueError("Normalized reference checksum mismatch")
        input_hashes[str(path)] = entry["normalized_sha256"]
        actual = subprocess.check_output(["bcftools", "query", "-l", str(path)], text=True).splitlines()
        if actual != samples:
            raise ValueError("Reference VCF sample order mismatch")
        region_file = reference_dir / chrom / "regions.tsv"
        if sha(region_file) != entry["regions_sha256"]:
            raise ValueError("Reference target-window checksum mismatch")
        input_hashes[str(region_file)] = entry["regions_sha256"]
        acquired = {tuple(line.split("\t")) for line in region_file.read_text().splitlines()}
        if any((chrom, row["start"], row["end"]) not in acquired for row in selected if row["chrom"] == chrom):
            raise ValueError("Reference acquisition does not cover the exact requested windows")
    cases, controls = gwas["case_count_from_portal"], gwas["control_count_from_portal"]
    nominal_n = cases + controls
    if nominal_n != gwas["total_participants_from_portal"]:
        raise ValueError("Inconsistent nominal external sample counts")
    for row in selected:
        input_hashes[str(source_dir/"regions"/f"{row['locus_id']}.tsv")] = source["region_files"][row["locus_id"]]["sha256"]
    out.mkdir(parents=True)

    def one(row):
        ident, chrom = row["locus_id"], row["chrom"]
        start, end = int(row["start"]), int(row["end"])
        folder = out / ident
        folder.mkdir()
        source_path = source_dir / "regions" / f"{ident}.tsv"
        if sha(source_path) != source["region_files"][ident]["sha256"]:
            raise ValueError(f"{ident}: source variant checksum mismatch")
        records = read_tsv(source_path)
        if any(r["locus_id"] != ident or not start <= int(r["pos"]) <= end for r in records):
            raise ValueError("Source variants lie outside their original window")
        normalized, norm_command = normalize_gwas(records, folder, chrom, fasta, contig_lengths[chrom])
        reference, query_command, ref_counts = read_reference(reference_dir / references[chrom]["normalized_file"], chrom, start, end, len(samples))
        matched, dosages, qc = match_variants(records, normalized, reference, start, end, nominal_n)
        write_tsv(folder / "variant_qc.tsv", qc, QC_FIELDS)
        write_tsv(folder / "matched_variants.tsv", matched, MATCH_FIELDS)
        exclusions = Counter(reason for r in qc if r["status"] == "excluded" for reason in r["reasons"].split(";"))
        by_class = {kind: dict(Counter(reason for r in qc if r["status"] == "excluded" and r["variant_type"] == kind
                                      for reason in r["reasons"].split(";"))) for kind in ["SNP", "INDEL", "unsupported"]}
        common = sum(common_af(r.get("af")) for r in records)
        status = "prepared" if len(matched) >= 2 else "fewer_than_two_matched_variants"
        info = dict(locus_id=ident, status=status, chrom=chrom, window_start=start, window_end=end,
            source_variants=len(records), source_common_variants=common, matched_variants=len(matched),
            matched_SNPs=sum(r["variant_type"] == "SNP" for r in matched),
            matched_INDELs=sum(r["variant_type"] == "INDEL" for r in matched),
            common_variant_coverage=len(matched)/common if common else None,
            exclusions_overlapping=dict(exclusions), exclusions_by_class_overlapping=by_class,
            reference_records_queried=sum(ref_counts.values()), reference_duplicate_keys=sum(v>1 for v in ref_counts.values()),
            normalization_command=norm_command, reference_query_command=query_command,
            gwas_input_sha256=sha(source_path), variant_qc_sha256=sha(folder/"variant_qc.tsv"),
            variants_file=str((folder/"matched_variants.tsv").relative_to(out)), variants_sha256=sha(folder/"matched_variants.tsv"),
            ld_variant_ids=[r["variant_id"] for r in matched])
        if len(matched) >= 2:
            matrix = signed_ld(dosages)
            ld_path = folder/"ld.tsv.gz"
            with gzip.open(ld_path, "wt") as stream:
                np.savetxt(stream, matrix, delimiter="\t", fmt="%.17g")
            info.update(ld_file=str(ld_path.relative_to(out)), ld_sha256=sha(ld_path), ld_shape=list(matrix.shape),
                        ld_symmetry_max_error=float(np.max(np.abs(matrix-matrix.T))),
                        ld_psd_basis="Gram correlation of complete observed ALT dosages; fit rechecks eigenvalues")
        (folder/"preparation.json").write_text(json.dumps(info, indent=2, allow_nan=False)+"\n")
        return info

    with ThreadPoolExecutor(max_workers=workers) as pool:
        results = list(pool.map(one, selected))
    if any(sha(path) != expected for path, expected in input_hashes.items()):
        raise ValueError("An input changed during preparation; outputs are incomplete")
    by_ident = {r["locus_id"]: r for r in results}
    qrows = []
    for row in leads:
        result = by_ident.get(row["locus_id"], {})
        qrows.append(dict(locus_id=row["locus_id"], chrom=row["chrom"], window_start=row["start"], window_end=row["end"],
            status=result.get("status", "not_selected_for_this_run"), source_variants=result.get("source_variants", "NA"),
            source_common_variants=result.get("source_common_variants", "NA"), matched_variants=result.get("matched_variants", "NA"),
            matched_SNPs=result.get("matched_SNPs", "NA"), matched_INDELs=result.get("matched_INDELs", "NA"),
            common_variant_coverage=result.get("common_variant_coverage", "NA"),
            exclusions_overlapping=json.dumps(result.get("exclusions_overlapping", {}), sort_keys=True)))
    write_tsv(out/"all_regions_qc.tsv", qrows, list(qrows[0]))
    manifest = dict(created_at_utc=datetime.now(timezone.utc).isoformat(), private_results=True,
        analysis="common_variant_external_LD_preparation_no_finemapping_fit", n_cases=cases, n_controls=controls,
        nominal_gwas_n=nominal_n, nominal_n_note="Trait-level portal count; variant-specific N is unavailable",
        provenance=dict(gwas_source=gwas.get("download_url", "External ALT-oriented GWAS summary statistics"),
            ld_source=ref.get("release", "Original unrelated EUR reference panel"), genome_build="GRCh38", effect_allele="ALT",
            ld_dosage_allele="ALT", allele_alignment="exact_normalized_REF_ALT", reference_missingness="complete_calls_only"),
        gwas_source_metadata=gwas, ld_source_metadata=ref, input_sha256=input_hashes, inputs_stable_during_read=True,
        source_dir=str(source_dir), reference_dir=str(reference_dir), source_manifest_sha256=sha(source_dir/"manifest.json"),
        reference_manifest_sha256=sha(reference_dir/"manifest.json"), fasta_sha256=ref["fasta_sha256"],
        script_sha256=sha(__file__), numpy_version=np.__version__,
        bcftools_version=subprocess.check_output(["bcftools", "--version"], text=True).splitlines()[:2],
        parameters=dict(min_GWAS_MAF=.01, min_reference_MAF=.01, max_absolute_AF_difference=.1,
            workers=workers, selected_locus=locus, p_value_filter=False, imputation=False,
            ranked_variant_truncation=False, normalized_duplicate_policy="exclude_all_ambiguous_keys"),
        source_regions=len(leads), requested_regions=len(selected), prepared_regions=sum(r["status"] == "prepared" for r in results),
        loci=[r for r in results if r["status"] == "prepared"],
        excluded_loci=[r for r in results if r["status"] != "prepared"],
        all_regions_qc_sha256=sha(out/"all_regions_qc.tsv"),
        limitations=["Common-variant preparation limited to the original lead-centered windows; no claim of whole-locus coverage.",
            "External European reference LD is an approximation to UKB LD and requires LD/z and model diagnostics before interpretation.",
            "All matching variants are retained regardless of association p-value; no variants are added by imputation.",
            "A missing haplotype/cohort association estimate is not repaired by this separate external-GWAS analysis.",
            "Exclusion counts overlap when a source variant fails multiple observable criteria."])
    (out/"input_manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False)+"\n")
    (out/"COMPLETED.txt").write_text("External LD preparation completed. No cohort genotype calls changed. No model fit by this step.\n")
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["source-dir", "reference-dir", "fasta", "out"]:
        parser.add_argument("--"+name, type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--locus")
    args = parser.parse_args()
    result = prepare(args.source_dir, args.reference_dir, args.fasta, args.out, args.workers, args.locus)
    print(json.dumps({key: result[key] for key in ["source_regions", "requested_regions", "prepared_regions"]}, indent=2))


if __name__ == "__main__":
    main()
