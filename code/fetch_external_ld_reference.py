#!/usr/bin/env python3
"""Acquire targeted, public GRCh38 reference genotypes for external GWAS LD.

This never reads cohort genotypes. The named IGSR phased release is already
filtered upstream and encodes retained FILTER as '.'. That release-specific
convention is recorded explicitly and does not change cohort PASS requirements.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import subprocess
import urllib.request

from extract_external_gwas import fingerprint, read_windows

BASE = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/"
PANEL = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def verify_input_snapshots(snapshots):
    """Require both content and metadata to match the inputs read for this run."""
    for name, before in snapshots.items():
        if fingerprint(before["path"]) != before:
            raise ValueError(f"Input changed during reference acquisition: {name}; output is incomplete")


def acquire(leads_path, fasta, out, workers=4):
    leads_path, fasta, out = (Path(p).resolve() for p in [leads_path, fasta, out])
    if out.exists():
        raise ValueError("Choose a new output directory")
    if not isinstance(workers, int) or not 1 <= workers <= 8:
        raise ValueError("Use 1..8 download workers")
    fasta_index = Path(str(fasta) + ".fai")
    if not fasta.is_file() or not fasta_index.is_file():
        raise ValueError("An indexed GRCh38 FASTA is required")
    snapshots = {name: fingerprint(path) for name, path in
                 [("leads", leads_path), ("fasta", fasta), ("fasta_index", fasta_index), ("script", Path(__file__))]}
    leads, _windows, _starts = read_windows(leads_path)
    out.mkdir(parents=True)
    panel = urllib.request.urlopen(PANEL, timeout=60).read()
    (out / "source_sample_panel.tsv").write_bytes(panel)
    samples = [r["sample"] for r in csv.DictReader(io.StringIO(panel.decode()), delimiter="\t") if r["super_pop"] == "EUR"]
    if len(samples) != 503 or len(set(samples)) != 503:
        raise ValueError("The original unrelated EUR reference panel must have 503 unique samples")
    sample_path = out / "EUR_503.samples.txt"
    sample_path.write_text("\n".join(samples) + "\n")

    def fetch(chrom):
        folder = out / chrom
        folder.mkdir()
        regions = folder / "regions.tsv"
        regions.write_text("".join(f"{r['chrom']}\t{r['start']}\t{r['end']}\n" for r in leads if r["chrom"] == chrom))
        url = BASE + f"1kGP_high_coverage_Illumina.{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
        raw, normalized = folder / "source_subset.vcf.gz", folder / "normalized.vcf.gz"
        commands = [
            ["bcftools", "view", "--no-version", "--regions-overlap", "0", "-R", str(regions), "-S", str(sample_path),
             "-m2", "-M2", "-v", "snps,indels", "-f", "PASS,.", "-Oz", "-o", str(raw), url],
            ["bcftools", "norm", "--no-version", "--check-ref", "e", "-f", str(fasta), "-Oz", "-o", str(normalized), str(raw)],
            ["bcftools", "index", str(normalized)],
        ]
        with (folder / "commands.log").open("w") as log:
            for command in commands:
                subprocess.run(command, cwd=folder, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=1800)
        actual = subprocess.check_output(["bcftools", "query", "-l", str(normalized)], text=True).splitlines()
        if actual != samples:
            raise ValueError(f"{chrom}: reference sample order changed")
        result = {"chrom": chrom, "source_url": url, "commands": commands,
                  "sample_count": len(actual), "source_subset_sha256": sha(raw),
                  "normalized_file": str(normalized.relative_to(out)), "normalized_sha256": sha(normalized),
                  "regions_sha256": sha(regions)}
        (folder / "manifest.json").write_text(json.dumps(result, indent=2) + "\n")
        print(json.dumps({"chrom": chrom, "status": "complete", "bytes": normalized.stat().st_size}), flush=True)
        return result

    chromosomes = sorted({r["chrom"] for r in leads}, key=lambda c: int(c[3:]))
    with ThreadPoolExecutor(max_workers=workers) as pool:
        results = list(pool.map(fetch, chromosomes))
    verify_input_snapshots(snapshots)
    manifest = {"created_at_utc": datetime.now(timezone.utc).isoformat(), "public_reference_data": True,
                "purpose": "external_GWAS_LD_not_cohort_genotype_recovery", "genome_build": "GRCh38",
                "release": BASE, "sample_panel_url": PANEL, "sample_panel_sha256": sha(out / "source_sample_panel.tsv"),
                "sample_count": 503, "population": "EUR_original_2504_unrelated_subset",
                "upstream_filtered_release_ACCEPTS_PASS_or_dot": True,
                "normalization": "bcftools norm; reference mismatches are fatal; no REF/ALT swaps",
                "fasta_sha256": snapshots["fasta"]["sha256"], "lead_table_sha256": snapshots["leads"]["sha256"],
                "input_snapshots": snapshots, "inputs_unchanged_during_acquisition": True,
                "script_sha256": snapshots["script"]["sha256"],
                "bcftools_version": subprocess.check_output(["bcftools", "--version"], text=True).splitlines()[:2],
                "chromosomes": results}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (out / "COMPLETED.txt").write_text("Targeted external reference acquisition completed. No cohort calls changed.\n")
    return manifest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--leads", type=Path, required=True)
    p.add_argument("--fasta", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--workers", type=int, default=4)
    args = p.parse_args()
    acquire(args.leads, args.fasta, args.out, args.workers)


if __name__ == "__main__":
    main()
