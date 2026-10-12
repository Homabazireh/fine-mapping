#!/usr/bin/env python3
"""Compare SNP/indel haplotypes around GWAS leads across sample groups.

Inputs: a multisample VCF, a GWAS TSV, and an exact sample/group TSV.
Start at run() to follow the five analysis steps. The VCF is read once and
may be plain text, ordinary gzip, or BGZF.

TSV coordinates are 1-based; BED output is 0-based half-open. Requires Python
>=3.10, with SciPy only for --test-associations and matplotlib only for --plot.
"""
from __future__ import annotations

import argparse
import bisect
from collections import Counter, defaultdict
import csv
from dataclasses import dataclass, field
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
import gzip
import hashlib
from itertools import combinations
import json
import math
from pathlib import Path
import re
import sys
from typing import TextIO

from project_config import ConfigArgumentParser, DEFAULT_PATHS, resolve_project_path

VERSION = "1.0.1"
DEFAULT_GWAS = resolve_project_path(DEFAULT_PATHS['gwas'])
DEFAULT_VCF = resolve_project_path(DEFAULT_PATHS['vcf'])
DNA = set("ACGT")


def open_text(path: str | Path) -> TextIO:
    """Detect compression by magic bytes, not by the filename alone."""
    with open(path, "rb") as handle:
        magic = handle.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return open(path, encoding="utf-8-sig", newline="")


def writer(path: Path, columns: list[str]):
    handle = (gzip.open(path, "wt", newline="", compresslevel=1)
              if path.suffix == ".gz" else open(path, "w", newline=""))
    obj = csv.DictWriter(handle, columns, delimiter="\t", lineterminator="\n")
    obj.writeheader()
    return handle, obj


def write_table(path: Path, rows: list[dict], columns: list[str]):
    handle, obj = writer(path, columns)
    with handle:
        obj.writerows(rows)


def value_or_na(value):
    if value is None:
        return "NA"
    if isinstance(value, float):
        return "NA" if math.isnan(value) else f"{value:.8g}"
    return value


def fraction(numerator: int, denominator: int):
    return numerator / denominator if denominator else None


def chrom_number(chrom: str) -> int | None:
    text = re.sub(r"^chr", "", chrom, flags=re.I)
    return int(text) if text.isdigit() and 1 <= int(text) <= 22 else None


def vcf_header(path: str | Path):
    """Read sample IDs, contig lengths, and genome-build declarations."""
    samples, contigs, reference_lines = [], {}, []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("##contig=<"):
                ident = re.search(r"(?:<|,)ID=([^,>]+)", line)
                length = re.search(r"(?:<|,)length=(\d+)", line, re.I)
                if ident:
                    contigs[ident[1]] = int(length[1]) if length else None
            if line.startswith("##reference=") or line.startswith("##assembly="):
                reference_lines.append(line.strip())
            if line.startswith("#CHROM\t"):
                samples = line.rstrip("\r\n").split("\t")[9:]
                break
            if not line.startswith("#"):
                break
    if not samples or len(samples) != len(set(samples)):
        raise ValueError("VCF needs a #CHROM header with unique sample IDs and genotype columns.")
    if not contigs:
        raise ValueError("VCF needs ##contig=<ID=...> header entries (length is optional).")
    return samples, contigs, reference_lines


def read_samples(path: str | Path, vcf_samples: list[str], allow_subset: bool = False):
    """No aliases, inferred phenotypes, or automatic group pooling."""
    mapping = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"sample", "group"} <= set(reader.fieldnames):
            raise ValueError("Sample TSV must have columns: sample and group.")
        for row_number, row in enumerate(reader, 2):
            sample, group = row.get("sample", ""), row.get("group", "")
            if sample is None or group is None or sample != sample.strip() or group != group.strip():
                raise ValueError(f"Whitespace or missing fields in sample table row {row_number}.")
            if not sample or group.upper() in {"", ".", "NA", "N/A", "UNRESOLVED", "UNKNOWN", "TODO"}:
                raise ValueError(f"Unresolved sample/group in sample table row {row_number}.")
            if sample in mapping:
                raise ValueError(f"Duplicate sample ID in metadata: {sample}")
            mapping[sample] = group
    extras, missing = set(mapping) - set(vcf_samples), set(vcf_samples) - set(mapping)
    if extras or (missing and not allow_subset):
        raise ValueError(f"Sample IDs do not match VCF exactly. Not in VCF: {sorted(extras)[:10]}; "
                         f"without group: {sorted(missing)[:10]}. "
                         "--allow-subset deliberately excludes unlisted VCF samples; it never renames them.")
    samples = [s for s in vcf_samples if s in mapping]
    groups = list(dict.fromkeys(mapping.values()))
    if len(groups) < 2 or not samples:
        raise ValueError("At least two explicitly defined nonempty groups are needed.")
    return samples, {s: mapping[s] for s in samples}, groups, sorted(missing)


@dataclass(frozen=True)
class Lead:
    chrom: str
    pos: int
    ea: str
    oa: str
    p: Decimal
    source_row: int
    z_abs: Decimal = Decimal(0)
    liftover_strand: str = "."
    input_id: str = "."

    def rank(self):
        return self.p, -self.z_abs, chrom_number(self.chrom), self.pos, self.ea, self.oa, self.source_row


ALIASES = {
    "chrom": ["Chromsome", "Chromosome", "chrom", "chr", "CHR", "CHROM", "#CHROM"],
    "pos": ["Position", "pos", "POS", "BP"],
    "p": ["Pval", "P", "p", "p_value", "pvalue", "P_VALUE"],
    "ea": ["EffectAllele", "effect_allele", "EA", "A1"],
    "oa": ["NonEffectAllele", "other_allele", "OA", "A2", "NEA"],
}


def read_gwas(args, contigs: dict):
    """Validate GWAS rows, map contig names, and retain P <= the cutoff."""
    rows, counts = [], Counter()
    threshold = Decimal(args.p_threshold)
    contigs_by_number = defaultdict(list)
    for contig in contigs:
        number = chrom_number(contig)
        if number is not None:
            contigs_by_number[number].append(contig)
    with open_text(args.gwas) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        columns = {}
        for key, aliases in ALIASES.items():
            explicit = getattr(args, f"{key}_col")
            column = explicit or next((name for name in aliases if name in fields), None)
            if column not in fields:
                raise ValueError(f"No GWAS {key} column. Use --{key}-col; available: {fields}")
            columns[key] = column
        seen = set()
        for row_number, row in enumerate(reader, 2):
            counts["input_rows"] += 1
            try:
                chrom = row[columns["chrom"]].strip()
                pos = int(row[columns["pos"]])
                pv = Decimal(row[columns["p"]])
                ea, oa = row[columns["ea"]].upper().strip(), row[columns["oa"]].upper().strip()
            except (TypeError, ValueError, InvalidOperation, AttributeError) as exc:
                raise ValueError(f"Malformed GWAS row {row_number}: {exc}") from exc
            if not pv.is_finite() or not 0 <= pv <= 1 or pos < 1:
                raise ValueError(f"Invalid position or p-value at GWAS row {row_number}.")
            if pv == 0:
                counts["zero_p_values_cannot_be_ranked_by_unavailable_precision"] += 1
            number = chrom_number(chrom)
            if number is None:
                counts["non_autosomal_rows_skipped"] += 1
                continue
            if chrom not in contigs:
                alternatives = contigs_by_number[number]
                if len(alternatives) != 1:
                    raise ValueError(f"Cannot map GWAS contig {chrom!r} uniquely to VCF header.")
                chrom = alternatives[0]
                counts["rows_with_contig_name_translation"] += 1
            if contigs[chrom] and pos > contigs[chrom]:
                raise ValueError(f"GWAS position outside contig at row {row_number}.")
            if pv > threshold:
                counts["above_p_threshold"] += 1
                continue
            if not ea or not oa or ea == oa or (set(ea) | set(oa)) - DNA:
                counts["non_sequence_or_invalid_allele_rows_skipped"] += 1
                continue
            key = chrom, pos, ea, oa
            if key in seen:
                raise ValueError(f"Duplicate significant GWAS allele key at row {row_number}: {key}")
            seen.add(key)
            z_abs = Decimal(0)
            if row.get("Beta") not in (None, "", ".", "NA") and row.get("SE") not in (None, "", ".", "NA"):
                try:
                    beta, se = Decimal(row["Beta"]), Decimal(row["SE"])
                    if beta.is_finite() and se.is_finite() and se > 0:
                        z_abs = abs(beta / se)
                except InvalidOperation:
                    pass  # Optional tie-breaker only; p-values remain validated.
            strand = row.get("liftover_strand", ".") or "."
            if strand not in {".", "+", "-"}:
                raise ValueError(f"Invalid liftover_strand in GWAS row {row_number}.")
            rows.append(Lead(chrom, pos, ea, oa, pv, row_number, z_abs, strand,
                             row.get("SNP", row.get("rsid", row.get("ID", "."))) or "."))
    counts["eligible_significant_rows"] = len(rows)
    return rows, dict(counts), columns


def choose_leads(rows: list[Lead], gap: int):
    """Preserve the uploaded consecutive-significant-gap rule; this is NOT LD clumping."""
    clusters = []
    for row in sorted(rows, key=lambda r: (chrom_number(r.chrom), r.pos, r.source_row)):
        if not clusters or row.chrom != clusters[-1][-1].chrom or row.pos - clusters[-1][-1].pos > gap:
            clusters.append([])
        clusters[-1].append(row)
    return sorted([(min(cluster, key=lambda r: r.rank()), len(cluster), cluster[0].pos, cluster[-1].pos)
                   for cluster in clusters], key=lambda item: item[0].rank())


@dataclass(frozen=True)
class Call:
    gt: tuple[int, int] | None
    phased: bool = False
    ps: str | None = None
    reason: str = "ok"


def parse_call(text: str, ix: dict[str, int], n_alleles: int, min_dp: float, min_gq: float) -> Call:
    """Parse one diploid genotype and apply its depth/quality filters."""
    fields = text.split(":")

    def get(name):
        index = ix.get(name, -1)
        return fields[index] if 0 <= index < len(fields) else "."

    gt = get("GT")
    if "." in gt or not gt:
        return Call(None, reason="missing_or_partial_GT")
    parts = re.split(r"[/|]", gt)
    if len(parts) != 2:
        return Call(None, reason="not_diploid")
    if not all(part.isdigit() for part in parts):
        return Call(None, reason="invalid_GT")
    alleles = tuple(int(part) for part in parts)
    if any(allele >= n_alleles for allele in alleles):
        return Call(None, reason="invalid_allele_index")
    if get("FT") not in {".", "", "PASS"}:
        return Call(None, reason="genotype_filter_failed")
    for name, threshold in (("DP", min_dp), ("GQ", min_gq)):
        if threshold <= 0:  # An explicitly disabled filter does not require its field.
            continue
        try:
            value = float(get(name))
        except ValueError:
            return Call(None, reason=f"missing_{name}")
        if not math.isfinite(value) or value < threshold:
            return Call(None, reason=f"low_or_invalid_{name}")
    ps = get("PS")
    if ps in {"", "."}:
        ps = None
    return Call(alleles, "|" in gt, ps)


def classify_variant(ref: str, alts: tuple[str, ...], max_indel_length: int) -> str | None:
    """Sequence-resolved SNPs/small length-changing alleles; multiallelic sites stay intact."""
    if not ref or not alts or any(not a or set(a) - DNA for a in (ref,) + alts):
        return None  # Includes symbolic SVs, '*', breakends, <NON_REF>, and N alleles.
    types = set()
    for alt in alts:
        if alt == ref:
            return None
        if len(ref) == len(alt) == 1:
            types.add("SNP")
        elif len(ref) != len(alt) and max(len(ref), len(alt)) <= max_indel_length + 1:
            types.add("indel")
        else:
            return None  # MNVs, long alleles and other complex records are outside scope.
    return "SNP/indel" if len(types) > 1 else next(iter(types))


@dataclass
class Marker:
    chrom: str
    pos: int
    ref: str
    alts: tuple[str, ...]
    kind: str
    calls: list[Call]
    source_line: int

    @property
    def alleles(self):
        return (self.ref,) + self.alts

    @property
    def key(self):
        return f"{self.chrom}:{self.pos}:{self.ref}:{','.join(self.alts)}"

    @property
    def called(self):
        return sum(call.gt is not None for call in self.calls)

    @property
    def minor_count(self):
        # For multiallelic records this means all non-major alleles combined.
        counts = Counter(a for call in self.calls if call.gt is not None for a in call.gt)
        return sum(counts.values()) - max(counts.values()) if counts else 0


def resolve_diplotype(calls: list[Call]):
    """Return an unordered pair, invariant to flipping all phased GTs in a block.

    Homozygous sites need no PS. With <=1 heterozygous marker the unordered pair
    is uniquely determined by genotype, even without read-based phase. With >=2
    heterozygous markers every heterozygote must be phased with the same explicit
    PS. PS values are compared only within this participant and chromosome.
    """
    if not calls or any(call.gt is None for call in calls):
        return None, "missing_or_low_quality", None
    heterozygotes = [call for call in calls if call.gt[0] != call.gt[1]]
    if len(heterozygotes) > 1:
        if any(not call.phased for call in heterozygotes):
            return None, "unphased_heterozygote", None
        if any(call.ps is None for call in heterozygotes):
            return None, "missing_PS", None
        sets = {call.ps for call in heterozygotes}
        if len(sets) != 1:
            return None, "different_phase_sets", None
        ps, status = next(iter(sets)), "same_PS"
    else:
        ps = heterozygotes[0].ps if heterozygotes and heterozygotes[0].phased else None
        status = "single_het_resolved" if heterozygotes else "homozygous_resolved"
    pair = tuple(sorted((tuple(call.gt[0] for call in calls), tuple(call.gt[1] for call in calls))))
    return pair, status, ps


@dataclass
class Locus:
    ident: str
    lead: Lead
    start: int
    end: int
    cluster_rows: int
    first_sig: int
    last_sig: int
    nearby: list[Marker] = field(default_factory=list)
    anchors: list[Marker] = field(default_factory=list)
    counts: Counter = field(default_factory=Counter)
    status: str = "pending"


def matching_anchor(marker: Marker, lead: Lead) -> bool:
    return marker.pos == lead.pos and lead.ea in marker.alleles and lead.oa in marker.alleles


def is_palindromic(lead: Lead):
    return len(lead.ea) == len(lead.oa) == 1 and {lead.ea, lead.oa} in ({"A", "T"}, {"C", "G"})


class WindowIndex:
    def __init__(self, loci: list[Locus]):
        self.entries = defaultdict(list)
        self.starts, self.maxends = {}, {}
        for locus in loci:
            self.entries[locus.lead.chrom].append(locus)
        for chrom, entries in self.entries.items():
            entries.sort(key=lambda locus: (locus.start, locus.end))
            self.starts[chrom] = [locus.start for locus in entries]
            running, ends = 0, []
            for locus in entries:
                running = max(running, locus.end)
                ends.append(running)
            self.maxends[chrom] = ends

    def query(self, chrom: str, pos: int):
        result = []
        entries = self.entries.get(chrom, [])
        i = bisect.bisect_right(self.starts.get(chrom, []), pos) - 1
        while i >= 0 and self.maxends[chrom][i] >= pos:
            if entries[i].end >= pos:
                result.append(entries[i])
            i -= 1
        return result


REGIONAL_FIELDS = ["locus_id", "variant_id", "chrom", "pos", "ref", "alt", "type", "distance_to_lead_bp",
                   "is_exact_lead_match", "group", "n_group", "n_called", "call_rate", "allele_index",
                   "allele", "allele_copies", "called_chromosomes", "allele_frequency", "carriers", "carrier_frequency"]
QC_FIELDS = ["locus_id", "variant_id", "type", "n_called", "n_samples", "call_rate", "nonmajor_allele_count",
             "eligible_nearby_marker", "call_exclusion_reasons"]


def scan_vcf(args, vcf_samples, samples, mapping, groups, loci, out):
    """Read once; retain only each locus's lead and nearest K eligible markers in memory."""
    index = WindowIndex(loci)
    sample_indices = [vcf_samples.index(sample) for sample in samples]
    group_indices = {group: [i for i, sample in enumerate(samples) if mapping[sample] == group] for group in groups}
    counts, seen_keys, visited = Counter(), set(), set()
    last_chrom, last_pos = None, 0
    cached_format, ix = None, {}
    rh, rw = writer(out / "regional_variant_frequencies.tsv.gz", REGIONAL_FIELDS)
    qh, qw = writer(out / "regional_variant_qc.tsv.gz", QC_FIELDS)
    with rh, qh, open_text(args.vcf) as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith("#") or not line.strip():
                continue
            # Avoid parsing sample columns for variants outside all selected windows.
            fields = line.rstrip("\r\n").split("\t", 9)
            if len(fields) != 10:
                raise ValueError(f"Malformed VCF line {line_number}: expected genotype columns.")
            chrom, pos = fields[0], int(fields[1])
            if pos < 1:
                raise ValueError(f"Invalid VCF position at line {line_number}.")
            if chrom != last_chrom:
                if chrom in visited:
                    raise ValueError(f"VCF is not coordinate sorted (contig reappears): {chrom}")
                visited.add(chrom)
                last_chrom, last_pos = chrom, 0
            if pos < last_pos:
                raise ValueError(f"VCF is not coordinate sorted at line {line_number}.")
            last_pos = pos
            counts["vcf_records_scanned"] += 1
            hits = index.query(chrom, pos)
            if not hits:
                continue
            counts["records_in_windows"] += 1
            for locus in hits:
                locus.counts["records_in_window"] += 1
                if pos == locus.lead.pos:
                    locus.counts["records_at_lead_position"] += 1
            ref, alts = fields[3], tuple(fields[4].split(","))
            key = chrom, pos, ref, alts
            if key in seen_keys:
                raise ValueError(f"Duplicate VCF variant record in selected windows at {key}; deduplicate upstream.")
            seen_keys.add(key)
            kind = classify_variant(ref, alts, args.max_indel_length)
            if kind is None:
                counts["excluded_not_SNP_or_small_indel"] += 1
                continue
            if fields[6] != "PASS" and not (args.allow_filter_dot and fields[6] == "."):
                counts["excluded_site_FILTER"] += 1
                continue
            cells = fields[9].split("\t")
            if len(cells) != len(vcf_samples):
                raise ValueError(f"VCF sample count mismatch at line {line_number}.")
            if cached_format != fields[8]:
                cached_format = fields[8]
                names = cached_format.split(":")
                if len(names) != len(set(names)) or "GT" not in names:
                    raise ValueError(f"Invalid FORMAT at VCF line {line_number}.")
                ix = {name: i for i, name in enumerate(names)}
            calls = [parse_call(cells[i], ix, len(alts) + 1, args.min_dp, args.min_gq) for i in sample_indices]
            if any(call.reason in {"invalid_GT", "invalid_allele_index"} for call in calls):
                raise ValueError(f"Invalid genotype at VCF line {line_number}.")
            marker = Marker(chrom, pos, ref, alts, kind, calls, line_number)
            call_rate = marker.called / len(samples)
            eligible = call_rate >= args.min_call_rate and marker.minor_count >= args.min_mac
            counts["passing_SNP_indel_records"] += 1
            counts["high_quality_homozygous_reference_genotypes"] += sum(call.gt == (0, 0) for call in calls)
            counts["high_quality_genotypes"] += marker.called
            by_group = {}
            for group, indices in group_indices.items():
                gt = [calls[i].gt for i in indices if calls[i].gt is not None]
                by_group[group] = (gt, Counter(a for g in gt for a in g))
            for locus in hits:
                locus.counts["passing_SNP_indel_records"] += 1
                exact = matching_anchor(marker, locus.lead)
                if exact:
                    locus.anchors.append(marker)
                # Other records at the anchor coordinate are not used as substitute leads/flanks.
                if eligible and pos != locus.lead.pos:
                    locus.nearby.append(marker)
                    locus.nearby.sort(key=lambda m: (abs(m.pos - locus.lead.pos), m.pos, m.key))
                    del locus.nearby[args.max_markers - 1:]
                qw.writerow(dict(locus_id=locus.ident, variant_id=marker.key, type=kind, n_called=marker.called,
                                 n_samples=len(samples), call_rate=value_or_na(call_rate),
                                 nonmajor_allele_count=marker.minor_count, eligible_nearby_marker=int(eligible),
                                 call_exclusion_reasons=json.dumps(dict(Counter(c.reason for c in calls if c.gt is None)), sort_keys=True)))
                for group in groups:
                    gt, allele_counts = by_group[group]
                    for allele_index, allele in enumerate(marker.alleles):
                        carriers = sum(allele_index in g for g in gt)
                        rw.writerow(dict(locus_id=locus.ident, variant_id=marker.key, chrom=chrom, pos=pos,
                                         ref=ref, alt=",".join(alts), type=kind, distance_to_lead_bp=pos - locus.lead.pos,
                                         is_exact_lead_match=int(exact), group=group, n_group=len(group_indices[group]),
                                         n_called=len(gt), call_rate=value_or_na(fraction(len(gt), len(group_indices[group]))),
                                         allele_index=allele_index, allele=allele, allele_copies=allele_counts[allele_index],
                                         called_chromosomes=2 * len(gt), allele_frequency=value_or_na(fraction(allele_counts[allele_index], 2 * len(gt))),
                                         carriers=carriers, carrier_frequency=value_or_na(fraction(carriers, len(gt)))))
    return dict(counts)


def bh_adjust(p_values: list[float | None]):
    """One Benjamini-Hochberg family across every finite haplotype carrier test."""
    valid = sorted((p, i) for i, p in enumerate(p_values) if p is not None and math.isfinite(p))
    result, running = [None] * len(p_values), 1.0
    total = len(valid)
    for rank in range(total, 0, -1):
        p, index = valid[rank - 1]
        running = min(running, p * total / rank)
        result[index] = running
    return result


def carrier_test(a: int, n_a: int, b: int, n_b: int):
    """Participant-level 2x2 test, never treating two chromosomes as independent people."""
    from scipy.stats import fisher_exact
    from scipy.stats.contingency import odds_ratio
    table = [[a, n_a - a], [b, n_b - b]]
    p = float(fisher_exact(table, alternative="two-sided").pvalue)
    result = odds_ratio(table, kind="conditional")
    interval = result.confidence_interval(confidence_level=0.95)
    return float(result.statistic), float(interval.low), float(interval.high), p


MARKER_FIELDS = ["locus_id", "marker_number", "variant_id", "chrom", "pos", "ref", "alt", "type", "is_lead",
                 "allele_index_to_sequence", "lead_effect_allele_index", "pooled_call_rate", "nonmajor_allele_count"]
SAMPLE_FIELDS = ["locus_id", "sample", "group", "status", "PS", "n_markers", "haplotype_1", "haplotype_2",
                 "haplotype_1_alleles", "haplotype_2_alleles"]
FREQ_FIELDS = ["locus_id", "haplotype_id", "allele_indices", "allele_sequences", "carries_lead_effect_allele", "group",
               "n_group", "n_resolved", "n_unresolved", "resolved_fraction", "haplotype_copies", "resolved_chromosomes",
               "haplotype_frequency", "carriers", "carrier_frequency"]
COMP_FIELDS = ["locus_id", "haplotype_id", "group_a", "group_b", "n_resolved_a", "n_resolved_b", "carriers_a", "carriers_b",
               "carrier_frequency_a", "carrier_frequency_b", "test_status", "conditional_odds_ratio_a_vs_b",
               "CI95_low", "CI95_high", "p_fisher_two_sided", "q_BH_all_haplotype_tests"]
GROUP_QC_FIELDS = ["locus_id", "group", "n_group", "n_resolved", "n_unresolved", "resolved_fraction", "exclusion_reasons"]
GENOTYPE_FIELDS = ["locus_id", "sample", "group", "variant_id", "GT", "PS", "call_status"]


def analyse_loci(args, samples, mapping, groups, loci, out):
    """Resolve each shared marker panel and write group counts and comparisons."""
    marker_rows, sample_rows, freq_rows, comparison_rows, group_qc, genotype_rows = [], [], [], [], [], []
    plot_payload = []
    group_sizes = Counter(mapping.values())
    for locus in loci:
        if locus.status != "pending":
            continue
        if not locus.anchors:
            locus.status = "no_unique_PASS_sequence_allele_match_at_lead"
            continue
        if len(locus.anchors) > 1:
            locus.status = "multiple_matching_lead_records"
            continue
        anchor = locus.anchors[0]
        if anchor.called / len(samples) < args.min_call_rate:
            locus.status = "lead_below_min_call_rate"
            continue
        panel = sorted([anchor] + locus.nearby, key=lambda marker: (marker.pos, marker.key))
        if len(panel) < 2:
            locus.status = "fewer_than_two_eligible_markers"
            continue
        locus.status = "analysed"
        locus.counts["selected_markers"] = len(panel)
        anchor_index, ea_index = panel.index(anchor), anchor.alleles.index(locus.lead.ea)
        for number, marker in enumerate(panel, 1):
            marker_rows.append(dict(locus_id=locus.ident, marker_number=number, variant_id=marker.key,
                                    chrom=marker.chrom, pos=marker.pos, ref=marker.ref, alt=",".join(marker.alts), type=marker.kind,
                                    is_lead=int(marker is anchor), allele_index_to_sequence=json.dumps(dict(enumerate(marker.alleles))),
                                    lead_effect_allele_index=ea_index if marker is anchor else "NA",
                                    pooled_call_rate=value_or_na(marker.called / len(samples)), nonmajor_allele_count=marker.minor_count))
        assignments, pooled = {}, Counter()
        for i, sample in enumerate(samples):
            calls = [marker.calls[i] for marker in panel]
            pair, status, ps = resolve_diplotype(calls)
            assignments[sample] = pair, status, ps
            if pair:
                pooled.update(pair)
            for marker, call in zip(panel, calls):
                gt = (("|" if call.phased else "/").join(map(str, call.gt)) if call.gt else ".")
                genotype_rows.append(dict(locus_id=locus.ident, sample=sample, group=mapping[sample], variant_id=marker.key,
                                          GT=gt, PS=call.ps if call.ps is not None else ".", call_status=call.reason))
        # Phenotype-blind deterministic IDs for identical allele strings at the SAME markers.
        ids = {hap: f"{locus.ident}_H{i:03d}" for i, hap in enumerate(sorted(pooled, key=lambda h: (-pooled[h], h)), 1)}
        sequence = lambda hap: "|".join(marker.alleles[a] for marker, a in zip(panel, hap))
        for sample, (pair, status, ps) in assignments.items():
            sample_rows.append(dict(locus_id=locus.ident, sample=sample, group=mapping[sample], status=status,
                                    PS=ps if ps is not None else ".", n_markers=len(panel),
                                    haplotype_1=ids[pair[0]] if pair else "NA", haplotype_2=ids[pair[1]] if pair else "NA",
                                    haplotype_1_alleles=sequence(pair[0]) if pair else "NA",
                                    haplotype_2_alleles=sequence(pair[1]) if pair else "NA"))
        resolved = {}
        for group in groups:
            selected = [s for s in samples if mapping[s] == group]
            resolved[group] = [assignments[s][0] for s in selected if assignments[s][0] is not None]
            n = len(resolved[group])
            reasons = Counter(assignments[s][1] for s in selected if assignments[s][0] is None)
            group_qc.append(dict(locus_id=locus.ident, group=group, n_group=group_sizes[group], n_resolved=n,
                                 n_unresolved=group_sizes[group] - n, resolved_fraction=value_or_na(n / group_sizes[group]),
                                 exclusion_reasons=json.dumps(dict(reasons), sort_keys=True)))
        for hap in ids:
            carrier_counts = {}
            for group in groups:
                pairs = resolved[group]
                n, copies = len(pairs), sum(pair.count(hap) for pair in pairs)
                carriers = sum(hap in pair for pair in pairs)
                carrier_counts[group] = carriers
                freq_rows.append(dict(locus_id=locus.ident, haplotype_id=ids[hap], allele_indices="|".join(map(str, hap)),
                                      allele_sequences=sequence(hap), carries_lead_effect_allele=int(hap[anchor_index] == ea_index), group=group,
                                      n_group=group_sizes[group], n_resolved=n, n_unresolved=group_sizes[group] - n,
                                      resolved_fraction=value_or_na(n / group_sizes[group]), haplotype_copies=copies, resolved_chromosomes=2 * n,
                                      haplotype_frequency=value_or_na(fraction(copies, 2 * n)), carriers=carriers,
                                      carrier_frequency=value_or_na(fraction(carriers, n))))
            for group_a, group_b in combinations(groups, 2):
                n_a, n_b = len(resolved[group_a]), len(resolved[group_b])
                a, b = carrier_counts[group_a], carrier_counts[group_b]
                status, odds, low, high, p = "not_requested", None, None, None, None
                if args.test_associations:
                    status = "tested_unadjusted"
                    if min(n_a, n_b) < args.min_test_samples:
                        status = "too_few_resolved_samples"
                    elif min(n_a / group_sizes[group_a], n_b / group_sizes[group_b]) < args.min_resolved_fraction:
                        status = "low_resolved_fraction"
                    elif a + b == 0 or a + b == n_a + n_b:
                        status = "no_carrier_variation"
                    elif min(a + b, n_a + n_b - a - b) < args.min_carriers:
                        status = "too_few_carriers_or_noncarriers"
                    else:
                        odds, low, high, p = carrier_test(a, n_a, b, n_b)
                comparison_rows.append(dict(locus_id=locus.ident, haplotype_id=ids[hap], group_a=group_a, group_b=group_b,
                                            n_resolved_a=n_a, n_resolved_b=n_b, carriers_a=a, carriers_b=b,
                                            carrier_frequency_a=fraction(a, n_a), carrier_frequency_b=fraction(b, n_b), test_status=status,
                                            conditional_odds_ratio_a_vs_b=odds, CI95_low=low, CI95_high=high, p_fisher_two_sided=p))
        if args.plot and len(plot_payload) < args.max_plots:
            plot_payload.append((locus, panel, assignments))
    adjusted = bh_adjust([row["p_fisher_two_sided"] for row in comparison_rows])
    for row, q in zip(comparison_rows, adjusted):
        row["q_BH_all_haplotype_tests"] = q
        for key in COMP_FIELDS:
            row[key] = value_or_na(row[key])
    write_table(out / "selected_markers.tsv", marker_rows, MARKER_FIELDS)
    write_table(out / "sample_haplotypes.tsv", sample_rows, SAMPLE_FIELDS)
    write_table(out / "haplotype_group_frequencies.tsv", freq_rows, FREQ_FIELDS)
    write_table(out / "haplotype_comparisons.tsv", comparison_rows, COMP_FIELDS)
    write_table(out / "group_phase_qc.tsv", group_qc, GROUP_QC_FIELDS)
    write_table(out / "selected_marker_genotypes.tsv.gz", genotype_rows, GENOTYPE_FIELDS)
    for payload in plot_payload:
        plot_locus(*payload, samples, mapping, groups, out)
    return sum(row["test_status"] == "tested_unadjusted" for row in comparison_rows)


def plot_locus(locus, panel, assignments, samples, mapping, groups, out):
    """Paper-like marker matrix; two rows/person, same marker columns in all groups.

    Uses matplotlib's default colormap; plots allele INDEX, not risk/protection.
    Unresolved full-panel haplotypes are NA, never falsely stitched across PS.
    """
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ordered = [sample for group in groups for sample in samples if mapping[sample] == group]
    matrix, labels = [], []
    for sample in ordered:
        pair, _, _ = assignments[sample]
        matrix.extend(pair if pair else ([float("nan")] * len(panel), [float("nan")] * len(panel)))
        labels.append(f"{mapping[sample]} / {sample}")
    array = np.array(matrix, dtype=float)
    fig, ax = plt.subplots(figsize=(max(9, len(panel) * 0.9), max(5, len(ordered) * 0.19)))
    maximum = max(len(marker.alts) for marker in panel)
    image = ax.imshow(array, aspect="auto", interpolation="none", vmin=0, vmax=maximum)
    ax.set_xticks(range(len(panel)))
    ax.set_xticklabels([f"{'* ' if marker.pos == locus.lead.pos else ''}{marker.pos}\n{marker.ref}>{','.join(marker.alts)}"
                        for marker in panel], rotation=60, ha="right", fontsize=8)
    ax.set_yticks([2 * i + 0.5 for i in range(len(ordered))])
    ax.set_yticklabels(labels, fontsize=6)
    for i, sample in enumerate(ordered):
        if assignments[sample][0] is None:
            ax.text((len(panel) - 1) / 2, 2 * i + 0.5, "NA", ha="center", va="center", fontsize=5)
    cursor = 0
    for group in groups[:-1]:
        cursor += sum(mapping[s] == group for s in ordered)
        ax.axhline(2 * cursor - 0.5, linewidth=1)
    colorbar = fig.colorbar(image, ax=ax, shrink=0.45, ticks=range(maximum + 1))
    colorbar.set_label("VCF allele index: 0=REF; 1,2,...=ALT")
    ax.set_title(f"{locus.ident}: {locus.lead.chrom}:{locus.lead.pos} (* GWAS lead)\n"
                 "Two rows/person; canonical unordered copies, not maternal/paternal; NA=unresolved")
    ax.set_xlabel("Selected SNP/indel marker positions (equally spaced columns, not a genomic distance scale)")
    fig.tight_layout()
    folder = out / "plots"
    folder.mkdir(exist_ok=True)
    fig.savefig(folder / f"{locus.ident}_haplotypes.png", dpi=160)
    plt.close(fig)


def argument_parser():
    p = ConfigArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                             path_defaults={'vcf': 'vcf', 'gwas': 'gwas', 'samples': 'samples', 'out': 'out'},
                             analysis_defaults=True)
    p.add_argument("--vcf", default=str(DEFAULT_VCF),
                   help="Multisample .vcf or .vcf.gz; defaults to the configured raw VCF")
    p.add_argument("--gwas", default=str(DEFAULT_GWAS),
                   help="GWAS TSV with chromosome, position, P, EA and OA; defaults to the configured filtered GWAS")
    p.add_argument("--samples", help="Authoritative TSV with sample and group columns")
    p.add_argument("--out", help="New or empty output directory")
    p.add_argument("--genome-build", required=True, choices=["GRCh37", "GRCh38"],
                   help="Your explicit declaration that BOTH inputs already use this build; no liftover is performed")
    p.add_argument("--p-threshold", default="1e-8",
                   help="Keep GWAS rows with P <= this value, even if the input was already filtered")
    p.add_argument("--lead-gap", type=int, default=150000, help="Maximum gap between consecutive significant GWAS positions in a cluster")
    p.add_argument("--top-loci", type=int, default=20, help="Top clusters by lead P; 0=all clusters")
    p.add_argument("--flank", type=int, default=50000, help="Bases on EACH side of the top variant; marker POS must lie inside")
    p.add_argument("--max-markers", type=int, default=8, help="Lead plus nearest eligible SNP/indel markers in a common panel")
    p.add_argument("--min-dp", type=float, default=10, help="Minimum genotype DP; 0 explicitly disables this filter")
    p.add_argument("--min-gq", type=float, default=20, help="Minimum genotype GQ; 0 explicitly disables this filter")
    p.add_argument("--min-call-rate", type=float, default=0.90, help="Pooled high-quality call rate required for lead and panel markers")
    p.add_argument("--min-mac", type=int, default=2, help="Minimum non-major allele count for nearby panel markers; not imposed on lead")
    p.add_argument("--max-indel-length", type=int, default=49, help="Maximum small indel size; stored REF/ALT length <= this + 1")
    p.add_argument("--allow-filter-dot", action=argparse.BooleanOptionalAction, help="Explicitly also accept site FILTER='.' (default PASS only)")
    p.add_argument("--allow-subset", action=argparse.BooleanOptionalAction, help="Deliberately use only listed sample IDs; still reject unknown IDs")
    p.add_argument("--expected-per-group", type=int, help="Optional strict sample-count check, e.g. 49")
    p.add_argument("--gwas-alleles-forward-confirmed", action="store_true",
                   help="You have verified forward-reference GWAS allele orientation, allowing palindromic/negative-liftover leads")
    p.add_argument("--test-associations", action=argparse.BooleanOptionalAction, help="Optional exploratory unadjusted participant-level carrier Fisher tests")
    p.add_argument("--reference-genotypes-confirmed", action="store_true",
                   help="You have verified genuine reference/nonreference calls, not variant-only merges or missing-to-reference filling")
    p.add_argument("--min-test-samples", type=int, default=10)
    p.add_argument("--min-resolved-fraction", type=float, default=0.80)
    p.add_argument("--min-carriers", type=int, default=3, help="Minimum pooled carriers AND noncarriers in each pairwise test")
    p.add_argument("--plot", action=argparse.BooleanOptionalAction, help="Write actual selected-marker haplotype PNG matrices")
    p.add_argument("--max-plots", type=int, default=10)
    for key in ALIASES:
        p.add_argument(f"--{key}-col", help=f"Override automatic GWAS {key} column detection")
    return p


@dataclass
class CohortData:
    """Validated VCF header and exact sample-to-group assignments."""

    vcf_samples: list[str]
    contigs: dict[str, int | None]
    reference_lines: list[str]
    samples: list[str]
    sample_groups: dict[str, str]
    groups: list[str]
    excluded_samples: list[str]
    group_counts: Counter


@dataclass
class GwasSelection:
    """Selected lead clusters plus the counts needed for the run summary."""

    chosen: list[tuple[Lead, int, int, int]]
    cluster_count: int
    counts: dict
    columns: dict[str, str]


def validate_run_options(args):
    """Check thresholds and requirements before reading analysis inputs."""
    if args.max_markers < 2 or min(args.flank, args.lead_gap, args.top_loci, args.min_dp, args.min_gq, args.min_mac) < 0:
        raise ValueError("max-markers must be >=2; distances, quality thresholds, top-loci and min-mac must be nonnegative.")
    if not 0 < args.min_call_rate <= 1 or not 0 <= args.min_resolved_fraction <= 1:
        raise ValueError("Call-rate/resolved-fraction thresholds must be within [0,1], with call rate >0.")
    if not 1 <= args.max_indel_length <= 49 or min(args.min_test_samples, args.min_carriers, args.max_plots) < 1:
        raise ValueError("Small-indel length must be 1..49; test sample/carrier/plot limits must be positive.")
    if args.expected_per_group is not None and args.expected_per_group < 1:
        raise ValueError("expected-per-group must be positive.")
    threshold = Decimal(args.p_threshold)
    if not threshold.is_finite() or not 0 <= threshold <= 1:
        raise ValueError("p-threshold must be within [0,1].")
    if args.test_associations and not args.reference_genotypes_confirmed:
        raise ValueError("Association tests require --reference-genotypes-confirmed AFTER independent callability/genotyping checks. "
                         "Descriptive output needs no such flag; missing genotypes are never changed.")
    if args.test_associations:
        try:
            from scipy.stats.contingency import odds_ratio  # noqa: F401
        except ImportError as exc:
            raise ValueError("Install scipy>=1.10 for association tests.") from exc
    if args.plot:
        try:
            import matplotlib  # noqa: F401
        except ImportError as exc:
            raise ValueError("Install matplotlib for --plot.") from exc


def read_cohort(args) -> CohortData:
    """Validate the declared genome build and the supplied sample groups."""
    vcf_samples, contigs, reference_lines = vcf_header(args.vcf)
    declared = " ".join(reference_lines).lower()
    opposite = ("grch37", "hg19") if args.genome_build == "GRCh38" else ("grch38", "hg38")
    if any(label in declared for label in opposite):
        raise ValueError("VCF reference/assembly header conflicts with --genome-build.")
    samples, sample_groups, groups, excluded_samples = read_samples(args.samples, vcf_samples, args.allow_subset)
    group_counts = Counter(sample_groups.values())
    if args.expected_per_group is not None and any(n != args.expected_per_group for n in group_counts.values()):
        raise ValueError(f"Group counts are {dict(group_counts)}, not {args.expected_per_group} in each group.")
    return CohortData(vcf_samples, contigs, reference_lines, samples, sample_groups,
                      groups, excluded_samples, group_counts)


def select_gwas_leads(args, contigs: dict) -> GwasSelection:
    """Choose the lowest-P lead in each physical-distance cluster."""
    gwas_rows, gwas_counts, columns = read_gwas(args, contigs)
    print(f"GWAS: {Path(args.gwas).resolve()}\n"
          f"Eligible variants at P <= {args.p_threshold}: {len(gwas_rows):,} / {gwas_counts.get('input_rows', 0):,}\n"
          f"Above threshold: {gwas_counts.get('above_p_threshold', 0):,}")
    clusters = choose_leads(gwas_rows, args.lead_gap)
    chosen = clusters if args.top_loci == 0 else clusters[:args.top_loci]
    if not chosen:
        raise ValueError("No eligible significant autosomal GWAS leads found.")
    return GwasSelection(chosen, len(clusters), gwas_counts, columns)


def prepare_run(args, cohort: CohortData, selection: GwasSelection):
    """Create a fresh output directory and record samples and locus windows."""
    out = Path(args.out)
    if out.exists() and any(out.iterdir()):
        raise ValueError(f"Output directory is not empty: {out}. Use a new directory to avoid stale results.")
    out.mkdir(parents=True, exist_ok=True)
    before = {str(Path(path).resolve()): (Path(path).stat().st_size, Path(path).stat().st_mtime_ns)
              for path in (args.vcf, args.gwas, args.samples)}
    loci = []
    for number, (lead, cluster_size, first, last) in enumerate(selection.chosen, 1):
        end = min(lead.pos + args.flank, cohort.contigs[lead.chrom] or (lead.pos + args.flank))
        locus = Locus(f"L{number:04d}", lead, max(1, lead.pos - args.flank),
                      end, cluster_size, first, last)
        if not args.gwas_alleles_forward_confirmed:
            if lead.liftover_strand == "-":
                locus.status = "negative_liftover_allele_orientation_unconfirmed"
            elif is_palindromic(lead):
                locus.status = "palindromic_lead_orientation_unconfirmed"
        loci.append(locus)
    sample_rows = [{"sample": sample, "group": cohort.sample_groups[sample]}
                   for sample in cohort.samples]
    write_table(out / "samples_used.tsv", sample_rows, ["sample", "group"])
    write_regions_bed(out / "regions.bed", loci)
    return out, loci, before


def write_regions_bed(path: Path, loci: list[Locus]):
    """Merge overlapping windows for extraction; keep each analysis locus separate."""
    intervals = sorted((locus.lead.chrom, locus.start - 1, locus.end) for locus in loci)
    merged = []
    for chrom, start, end in intervals:
        if merged and chrom == merged[-1][0] and start <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], end)
        else:
            merged.append([chrom, start, end])
    with open(path, "w") as handle:
        for interval in merged:
            handle.write("\t".join(map(str, interval)) + "\n")


def write_lead_loci(out: Path, loci: list[Locus]):
    """Record each selected lead, including loci skipped during quality checks."""
    locus_rows = [dict(
        locus_id=locus.ident,
        chrom=locus.lead.chrom,
        lead_pos=locus.lead.pos,
        lead_effect_allele=locus.lead.ea,
        lead_other_allele=locus.lead.oa,
        lead_p=str(locus.lead.p),
        gwas_source_line=locus.lead.source_row,
        input_variant_id=locus.lead.input_id,
        start=locus.start,
        end=locus.end,
        cluster_significant_rows=locus.cluster_rows,
        cluster_first_significant_pos=locus.first_sig,
        cluster_last_significant_pos=locus.last_sig,
        status=locus.status,
        records_in_window=locus.counts["records_in_window"],
        records_at_lead_position=locus.counts["records_at_lead_position"],
        exact_passing_lead_matches=len(locus.anchors),
        selected_markers=locus.counts["selected_markers"],
    ) for locus in loci]
    write_table(out / "lead_loci.tsv", locus_rows, list(locus_rows[0]))


def interpretation_warnings(args, scan_counts: dict) -> list[str]:
    """State the limits that apply to these descriptive or exploratory results."""
    warnings = ["Descriptive frequencies are conditional on observed high-quality resolved calls, not guaranteed unbiased population estimates.",
                "Lead clustering uses physical distance, not LD; GWAS leads and observed haplotypes do not establish causality.",
                "Local copy 1/2 and H001 IDs are not maternal/paternal or globally interchangeable across loci.",
                "Sequence indels are compared by exact representation; both GWAS and VCF must already be harmonized/normalized.",
                "GWAS effect allele is NOT automatically the risk allele; effect direction is not interpreted.",
                "Genotype/phase missingness can bias comparisons; inspect group_phase_qc.tsv before interpretation."]
    if not args.reference_genotypes_confirmed:
        warnings.append("Reference-call provenance is not confirmed. No inferential association tests were requested/permitted.")
    if scan_counts.get("high_quality_homozygous_reference_genotypes", 0) == 0:
        warnings.append("No high-quality 0/0 calls in selected windows: investigate variant-only merging/reference-call absence.")
    if args.test_associations:
        warnings.append("Fisher carrier tests are exploratory, unadjusted for ancestry, age, sex, technical batch or relatedness. "
                        "They require independent participants; no progression, protective or causal claim is supported.")
    return warnings


def finish_run(args, cohort, selection, loci, out, before, scan_counts, n_tests):
    """Write provenance and mark completion only after checking input stability."""
    write_lead_loci(out, loci)
    stable = all((Path(path).stat().st_size, Path(path).stat().st_mtime_ns) == stats
                 for path, stats in before.items())
    warnings = interpretation_warnings(args, scan_counts)
    locus_status_counts = dict(Counter(locus.status for locus in loci))
    summary = dict(
        version=VERSION,
        completed_at_utc=datetime.now(timezone.utc).isoformat(),
        genome_build_user_declared=args.genome_build,
        reference_header=cohort.reference_lines,
        parameters=vars(args),
        group_counts=dict(cohort.group_counts),
        excluded_unlisted_vcf_samples=cohort.excluded_samples,
        gwas_column_mapping=selection.columns,
        gwas_counts=selection.counts,
        distance_clusters_available=selection.cluster_count,
        selected_loci=len(loci),
        locus_status_counts=locus_status_counts,
        scan_counts=scan_counts,
        haplotype_carrier_tests_performed=n_tests,
        inputs_stable_during_read=stable,
        warnings=warnings,
        input_files={path: {"bytes": stat[0], "mtime_ns": stat[1]} for path, stat in before.items()},
        small_input_sha256={str(Path(args.samples).resolve()): hashlib.sha256(Path(args.samples).read_bytes()).hexdigest()},
    )
    (out / "run_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (out / "INTERPRETATION.txt").write_text("\n".join(warnings) + "\n")
    if not stable:
        raise ValueError("An input changed during the run; discard outputs and rerun with stable inputs.")
    (out / "COMPLETED.txt").write_text("Completed. Read run_summary.json and lead_loci.tsv, including skipped loci.\n")
    print(f"Groups: {dict(cohort.group_counts)}\nSelected loci: {len(loci)} / {selection.cluster_count} distance clusters\n"
          f"Locus status: {locus_status_counts}\n"
          f"Carrier tests: {n_tests}\nOutputs: {out.resolve()}")
    if not any(l.status == "analysed" for l in loci):
        print("WARNING: No locus yielded a usable multi-marker panel. Inspect lead_loci.tsv; do not interpret this as no genetic effect.", file=sys.stderr)
    return summary


def run(args):
    """Validate inputs, select leads, prepare windows, compare calls, and save results."""
    # 1. Validate analysis options and exact sample/group assignments.
    validate_run_options(args)
    cohort = read_cohort(args)

    # 2. Rank significant GWAS leads within physical-distance clusters.
    selection = select_gwas_leads(args, cohort.contigs)

    # 3. Prepare a fresh output directory and the lead-centered windows.
    out, loci, before = prepare_run(args, cohort, selection)

    # 4. Read the VCF once, then compare resolved haplotypes by group.
    scan_counts = scan_vcf(args, cohort.vcf_samples, cohort.samples,
                           cohort.sample_groups, cohort.groups, loci, out)
    n_tests = analyse_loci(args, cohort.samples, cohort.sample_groups, cohort.groups, loci, out)

    # 5. Save locus status, input provenance, and the completion marker.
    return finish_run(args, cohort, selection, loci, out, before, scan_counts, n_tests)


def main(argv=None):
    args = argument_parser().parse_args(argv)
    try:
        run(args)
    except (ValueError, OSError, InvalidOperation, csv.Error) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
