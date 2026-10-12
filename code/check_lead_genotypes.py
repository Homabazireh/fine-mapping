#!/usr/bin/env python3
"""Check genotype completeness at GWAS leads using the VCF index.

This reads aggregate calls for every VCF sample. Optional metadata adds group
counts for exact or explicitly verified sample-ID matches. It does not assign haplotypes or replace
missing genotypes.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
from types import SimpleNamespace

from compare_haplotypes import (
    ALIASES, choose_leads, classify_variant, is_palindromic, parse_call,
    read_gwas, value_or_na, vcf_header, write_table,
)
from make_sample_template import _resolve_metadata_groups, read_metadata_groups
from project_config import DEFAULT_CONFIG, PROJECT_ROOT, load_config, resolve_project_path


DEFAULT_OUT = PROJECT_ROOT / 'data' / 'intermediate' / 'cohort_preflight'
PARAMETERS = ('p_threshold', 'lead_gap', 'top_loci', 'min_dp', 'min_gq',
              'min_call_rate', 'max_indel_length', 'allow_filter_dot')


def input_info(path, include_hash=False):
    """Record source size/time; hash small inputs, never the full cohort VCF."""
    path = Path(path).resolve()
    stat = path.stat()
    info = dict(path=str(path), size_bytes=stat.st_size, mtime_ns=stat.st_mtime_ns)
    if include_hash:
        info['sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
    return info


def record_counts(fields, n_samples, settings):
    """Count valid diploid calls before and after genotype quality filters."""
    if len(fields) != n_samples + 9:
        raise ValueError(f'VCF sample count mismatch at {fields[0]}:{fields[1]}')
    alleles = [fields[3], *fields[4].split(',')]
    indices = {name: index for index, name in enumerate(fields[8].split(':'))}
    # Raw calls ignore DP, GQ and genotype FT; missing/invalid GTs stay uncalled.
    raw_indices = {name: index for name, index in indices.items() if name != 'FT'}
    raw = [parse_call(cell, raw_indices, len(alleles), 0, 0) for cell in fields[9:]]
    quality = [parse_call(cell, indices, len(alleles), settings['min_dp'], settings['min_gq'])
               for cell in fields[9:]]
    called = sum(call.gt is not None for call in quality)
    return dict(
        raw_called_samples=sum(call.gt is not None for call in raw),
        high_quality_called_samples=called,
        high_quality_call_rate=called / n_samples if n_samples else None,
        raw_reference_genotypes=sum(call.gt == (0, 0) for call in raw),
        high_quality_reference_genotypes=sum(call.gt == (0, 0) for call in quality),
        high_quality_exclusion_reasons=dict(Counter(call.reason for call in quality if call.gt is None)),
    )


def matching_records(lead, records, settings):
    """Return exact allele matches, passing records, and a unique record to count."""
    matches = [fields for fields in records
               if {lead.ea, lead.oa} <= {fields[3], *fields[4].split(',')}]
    passing = [fields for fields in matches
               if (fields[6] == 'PASS' or (settings['allow_filter_dot'] and fields[6] == '.'))
               and classify_variant(fields[3], tuple(fields[4].split(',')), settings['max_indel_length'])]
    # Ambiguous matches have no single sample denominator; report counts as NA.
    selected = passing[0] if len(passing) == 1 else matches[0] if len(matches) == 1 else None
    return matches, passing, selected


def inspect_lead(number, lead, records, n_samples, settings):
    """Summarize an exact lead match without choosing a replacement variant."""
    matches, passing, selected = matching_records(lead, records, settings)
    counts = record_counts(selected, n_samples, settings) if selected is not None else {}
    call_rate = counts.get('high_quality_call_rate')
    passes_call_rate = len(passing) == 1 and call_rate is not None and call_rate >= settings['min_call_rate']
    if lead.liftover_strand == '-':
        status = 'negative_liftover_allele_orientation_unconfirmed'
    elif is_palindromic(lead):
        status = 'palindromic_lead_orientation_unconfirmed'
    elif not passing:
        status = 'no_unique_PASS_sequence_allele_match_at_lead'
    elif len(passing) > 1:
        status = 'multiple_matching_lead_records'
    elif not passes_call_rate:
        status = 'lead_below_min_call_rate'
    else:
        status = 'lead_passes_QC_panel_not_evaluated'
    return dict(
        locus_id=f'L{number:04d}', chrom=lead.chrom, lead_pos=lead.pos,
        effect_allele=lead.ea, other_allele=lead.oa, lead_p=str(lead.p),
        palindromic=int(is_palindromic(lead)), liftover_strand=lead.liftover_strand,
        records_at_position=len(records), exact_allele_matches=len(matches),
        site_PASS_matches=sum(fields[6] == 'PASS' for fields in matches),
        passing_sequence_matches=len(passing), site_filters=','.join(sorted({fields[6] for fields in matches})),
        n_samples=n_samples, raw_called_samples=counts.get('raw_called_samples'),
        high_quality_called_samples=counts.get('high_quality_called_samples'),
        high_quality_call_rate=call_rate,
        raw_reference_genotypes=counts.get('raw_reference_genotypes'),
        high_quality_reference_genotypes=counts.get('high_quality_reference_genotypes'),
        passes_min_call_rate=int(passes_call_rate), lead_status=status,
        high_quality_exclusion_reasons=json.dumps(counts.get('high_quality_exclusion_reasons', {}), sort_keys=True),
    )


def inspect_lead_groups(number, lead, records, group_indices, settings):
    """Count each source category using its resolved VCF sample indices."""
    _, _, selected = matching_records(lead, records, settings)
    rows = []
    for group, indices in group_indices.items():
        counts = {}
        if selected is not None:
            fields = selected[:9] + [selected[9 + index] for index in indices]
            counts = record_counts(fields, len(indices), settings)
        rows.append(dict(
            locus_id=f'L{number:04d}', chrom=lead.chrom, lead_pos=lead.pos,
            group=group, n_group=len(indices),
            high_quality_called_samples=counts.get('high_quality_called_samples'),
            high_quality_call_rate=counts.get('high_quality_call_rate'),
            raw_called_samples=counts.get('raw_called_samples'),
            raw_reference_genotypes=counts.get('raw_reference_genotypes'),
            high_quality_reference_genotypes=counts.get('high_quality_reference_genotypes'),
        ))
    return rows


def check_lead_genotypes(config=DEFAULT_CONFIG, out_dir=DEFAULT_OUT, metadata_path=None, id_map_path=None):
    """Query only selected lead positions and save aggregate QC with provenance."""
    if id_map_path is not None and metadata_path is None:
        raise ValueError('--id-map requires --metadata')
    executable = shutil.which('bcftools')
    if executable is None:
        raise RuntimeError('bcftools is required for indexed lead-site checks.')
    out = Path(out_dir)
    if out.exists() and any(out.iterdir()):
        raise ValueError(f'Output directory is not empty: {out}. Choose a new directory.')
    settings = load_config(config)
    parameters = {name: settings['analysis'][name] for name in PARAMETERS}
    inputs = {
        'vcf': input_info(settings['paths']['vcf']),
        'gwas': input_info(settings['paths']['gwas'], include_hash=True),
        'config': input_info(resolve_project_path(config), include_hash=True),
    }
    metadata = None
    if metadata_path is not None:
        inputs['metadata'] = input_info(metadata_path, include_hash=True)
        metadata = read_metadata_groups(metadata_path)
    if id_map_path is not None:
        inputs['id_map'] = input_info(id_map_path, include_hash=True)
    samples, contigs, reference_lines = vcf_header(inputs['vcf']['path'])
    group_indices = {}
    exact_counts, alias_counts = Counter(), Counter()
    if metadata is not None:
        resolved_groups = _resolve_metadata_groups(metadata, samples, id_map_path)
        group_indices = {group: [] for group in metadata.values()}
        for index, sample in enumerate(samples):
            if sample in resolved_groups:
                group = resolved_groups[sample]
                group_indices[group].append(index)
                (exact_counts if sample in metadata else alias_counts)[group] += 1
    gwas_args = SimpleNamespace(gwas=inputs['gwas']['path'], p_threshold=parameters['p_threshold'],
                               **{f'{name}_col': None for name in ALIASES})
    gwas, gwas_counts, columns = read_gwas(gwas_args, contigs)
    clusters = choose_leads(gwas, parameters['lead_gap'])
    chosen = clusters if parameters['top_loci'] == 0 else clusters[:parameters['top_loci']]
    if not chosen:
        raise ValueError('No eligible significant autosomal GWAS leads found.')
    regions = ','.join(f'{lead.chrom}:{lead.pos}-{lead.pos}' for lead, *_ in chosen)
    command = [executable, 'view', '--no-version', '-H', '-r', regions, inputs['vcf']['path']]
    query = subprocess.run(command, capture_output=True, text=True)
    if query.returncode:
        raise RuntimeError(f'Indexed VCF query failed: {query.stderr.strip()}')
    records = defaultdict(list)
    for line in query.stdout.splitlines():
        fields = line.split('\t')
        records[fields[0], int(fields[1])].append(fields)
    rows = [inspect_lead(number, lead, records[lead.chrom, lead.pos], len(samples), parameters)
            for number, (lead, *_) in enumerate(chosen, 1)]
    group_rows = [row for number, (lead, *_) in enumerate(chosen, 1)
                  for row in inspect_lead_groups(number, lead, records[lead.chrom, lead.pos],
                                                 group_indices, parameters)]
    stable = all((Path(info['path']).stat().st_size, Path(info['path']).stat().st_mtime_ns)
                 == (info['size_bytes'], info['mtime_ns']) for info in inputs.values())
    if not stable:
        raise ValueError('An input changed during preflight; rerun with stable inputs.')
    summary = dict(
        schema_version=1, scope='all_vcf_samples_no_phenotypes',
        created_at_utc=datetime.now(timezone.utc).isoformat(),
        n_samples=len(samples), n_leads=len(rows),
        leads_passing_call_rate=sum(row['passes_min_call_rate'] for row in rows),
        leads_without_exact_match=sum(row['exact_allele_matches'] == 0 for row in rows),
        max_high_quality_called=max(row['high_quality_called_samples'] or 0 for row in rows),
        max_high_quality_call_rate=max(row['high_quality_call_rate'] or 0 for row in rows),
        raw_reference_genotypes=sum(row['raw_reference_genotypes'] or 0 for row in rows),
        high_quality_reference_genotypes=sum(row['high_quality_reference_genotypes'] or 0 for row in rows),
        lead_status_counts=dict(Counter(row['lead_status'] for row in rows)),
        parameters=parameters, inputs=inputs, vcf_reference_headers=reference_lines,
        gwas_counts=gwas_counts, gwas_column_mapping=columns,
        bcftools_command=command, bcftools_stderr=query.stderr.strip(),
        inputs_stable_during_read=stable,
        notes=['All-sample totals count every VCF sample independently of phenotype labels.',
               'Only exact lead positions were queried. Nearby marker panels and haplotypes were not analysed.',
               'Raw counts require valid diploid GTs but ignore DP, GQ and genotype FT.',
               'Missing genotypes remain missing. No reference-genotype or allele-orientation attestation was made.'],
    )
    out.mkdir(parents=True, exist_ok=True)
    write_table(out / 'lead_genotype_qc.tsv',
                [{key: value_or_na(value) for key, value in row.items()} for row in rows], list(rows[0]))
    if metadata is not None:
        group_table = out / 'group_lead_genotype_qc.tsv'
        group_fields = ['locus_id', 'chrom', 'lead_pos', 'group', 'n_group',
                        'high_quality_called_samples', 'high_quality_call_rate', 'raw_called_samples',
                        'raw_reference_genotypes', 'high_quality_reference_genotypes']
        write_table(group_table, [{key: value_or_na(value) for key, value in row.items()}
                                 for row in group_rows], group_fields)
        matched_count = sum(len(indices) for indices in group_indices.values())
        summary['group_qc'] = dict(
            scope='exact_and_verified_id_matches' if id_map_path is not None else 'exact_metadata_id_matches_only',
            source_category_counts=dict(Counter(metadata.values())),
            exact_match_category_counts={group: exact_counts[group] for group in group_indices},
            verified_alias_category_counts={group: alias_counts[group] for group in group_indices},
            matched_category_counts={group: len(indices) for group, indices in group_indices.items()},
            exact_match_sample_count=sum(exact_counts.values()),
            verified_alias_sample_count=sum(alias_counts.values()),
            matched_sample_count=matched_count,
            unmatched_vcf_count=len(samples) - matched_count,
            unmatched_metadata_count=len(metadata) - matched_count,
            table=str(group_table.resolve()),
        )
        matching_note = ('Group QC uses exact IDs plus the supplied verified one-to-one alias map.'
                         if id_map_path is not None else 'Group QC uses exact metadata IDs only.')
        summary['notes'].append(matching_note + ' Unmatched IDs are not assigned to any category; this is not a subset haplotype analysis.')
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(f"Checked {summary['n_leads']} leads across {summary['n_samples']} samples. "
          f"Passing call rate: {summary['leads_passing_call_rate']}. Outputs: {out.resolve()}")
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', default=str(DEFAULT_CONFIG), help='Project settings JSON')
    parser.add_argument('--out', type=Path, default=DEFAULT_OUT, help='New or empty preflight output directory')
    parser.add_argument('--metadata', type=Path, help='Optional category TSV for exact-ID group QC; aliases are not inferred')
    parser.add_argument('--id-map', type=Path, help='Verified aliases TSV with metadata_sample and vcf_sample columns; requires --metadata')
    args = parser.parse_args(argv)
    try:
        check_lead_genotypes(args.config, args.out, metadata_path=args.metadata, id_map_path=args.id_map)
    except (OSError, ValueError, RuntimeError) as error:
        parser.error(str(error))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
