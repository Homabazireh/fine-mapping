#!/usr/bin/env python3
"""Build the one-figure aggregate publication from completed analysis outputs."""
from __future__ import annotations

import argparse
import csv
import hashlib
from html import escape
import json
from pathlib import Path
import shutil

from export_panel_c import export_panel_c


def read_table(path):
    with Path(path).open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def validate_publication_inputs(run_dir, evidence_dir):
    """Bind evidence to this stable, real cohort run before creating outputs."""
    run_dir, evidence_dir = Path(run_dir), Path(evidence_dir)
    summary = json.loads((run_dir / 'run_summary.json').read_text())
    if not (run_dir / 'COMPLETED.txt').is_file():
        raise ValueError('Only completed runs may be published')
    if summary.get('synthetic') is True or summary.get('inputs_stable_during_read') is not True:
        raise ValueError('Publication requires real data and stable inputs')
    audit = json.loads((evidence_dir / 'audit_summary.json').read_text())
    if audit.get('synthetic') is not False:
        raise ValueError('Evidence audit must explicitly identify real data')
    def check(path, expected):
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f'Publication provenance mismatch: {path.name}')
    for key, filename in [('cohort_summary', 'run_summary.json'), ('selected_leads', 'lead_loci.tsv')]:
        check(run_dir / filename, audit['provenance'][key]['sha256'])
    for filename in ['external_gwas_leads.tsv', 'legacy_locus_audit.tsv', 'candidate_variants.tsv']:
        check(evidence_dir / filename, audit['output_sha256'][filename])
    if audit['cohort_group_counts'] != summary['group_counts']:
        raise ValueError('Evidence and cohort group counts differ')
    return summary


def publish(run_dir, evidence_dir, out):
    run_dir, evidence_dir, out = map(Path, (run_dir, evidence_dir, out))
    summary = validate_publication_inputs(run_dir, evidence_dir)
    manifest = export_panel_c(run_dir, out)
    # Public provenance uses a run label, not a local account or filesystem path.
    manifest['source_run'] = run_dir.name
    (out / 'panel_c_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    aggregate_files = ['lead_loci.tsv', 'haplotype_group_frequencies.tsv',
                       'haplotype_overall_comparisons.tsv', 'haplotype_comparisons.tsv',
                       'group_phase_qc.tsv']
    for name in aggregate_files:
        shutil.copyfile(run_dir / name, out / name)
    evidence_files = ['external_gwas_leads.tsv', 'legacy_locus_audit.tsv', 'candidate_variants.tsv']
    for name in evidence_files:
        shutil.copyfile(evidence_dir / name, out / name)
    shutil.copyfile(evidence_dir / 'audit_summary.json', out / 'variant_evidence_audit.json')
    shutil.copyfile(evidence_dir / 'REPORT.txt', out / 'VARIANT_EVIDENCE.txt')
    evidence_files += ['variant_evidence_audit.json', 'VARIANT_EVIDENCE.txt']
    overall = read_table(run_dir / 'haplotype_overall_comparisons.tsv')
    pairs = read_table(run_dir / 'haplotype_comparisons.tsv')
    candidates = read_table(evidence_dir / 'candidate_variants.tsv')
    significant = [r for r in pairs if r['test_status'] == 'tested_unadjusted'
                   and r['q_BH_all_pairwise_tests'] != 'NA' and float(r['q_BH_all_pairwise_tests']) < .05]
    valid_panels = summary['locus_status_counts'].get('analysed', 0)
    overall_n = sum(r['test_status'] == 'tested_unadjusted' for r in overall)
    pair_n = sum(r['test_status'] == 'tested_unadjusted' for r in pairs)
    status_rows = []
    for lead in read_table(run_dir / 'lead_loci.tsv'):
        ident = lead['locus_id']
        tested = any(r['locus_id'] == ident and r['test_status'] == 'tested_unadjusted' for r in overall)
        status_rows.append({'locus_id': ident, 'chrom': lead['chrom'], 'lead_pos': lead['lead_pos'],
                            'panel_status': lead['status'],
                            'overall_comparison': 'tested' if tested else 'not_estimable_or_not_permitted',
                            'association_tests_requested': bool(summary['parameters']['test_associations']),
                            'reference_genotypes_confirmed': bool(summary['parameters']['reference_genotypes_confirmed'])})
    with (out / 'comparison_status.tsv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, list(status_rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(status_rows)
    counts = {'group_counts': summary['group_counts'], 'selected_regions': summary['selected_loci'],
              'usable_haplotype_panels': valid_panels, 'overall_tests_performed': overall_n,
              'pairwise_tests_performed': pair_n, 'significant_pairwise_tests_FDR_005': len(significant),
              'candidate_variants_supported_by_cohort_haplotypes': len(candidates),
              'comparison_available': bool(overall_n or pair_n),
              'locus_status_counts': summary['locus_status_counts'],
              'high_quality_homozygous_reference_calls': summary['scan_counts']['high_quality_homozygous_reference_genotypes']}
    (out / 'comparison_results.json').write_text(json.dumps(counts, indent=2) + '\n')
    if not valid_panels:
        finding = ('No haplotype can currently be identified as differing between Diabetes, Prediabetes and Control. '
                   'All selected regions failed the requirements for constructing a usable haplotype panel. '
                   'Haplotype frequencies, effect sizes, confidence intervals, overall and pairwise p-values, '
                   'and adjusted p-values are not estimable from this run. Header-only result tables mean '
                   'not tested, not a null result. The single panel C therefore retains its genotype-completeness '
                   'meaning; it is not relabeled as a haplotype-frequency figure.')
    else:
        finding = (f'{valid_panels} panels produced descriptive haplotype estimates; {overall_n} overall and '
                   f'{pair_n} pairwise tests were performed. {len(significant)} pairwise comparisons have BH q < 0.05. '
                   'See the tables for denominators, reasons for skipped tests, effect direction and uncertainty. '
                   'The single panel C measures genotype completeness.')
    parameters = summary['parameters']
    methods = ('The following methods are implemented for data that meet the eligibility criteria; '
               f'this run performed {overall_n} overall and {pair_n} pairwise tests. '
               'Use the same phenotype-blind marker panel within each locus. Report chromosome-copy frequencies '
               'descriptively. Test participant carrier status using the exact Fisher–Freeman–Halton test '
               'across all three groups and two-sided Fisher tests for each pair. Report conditional odds ratios '
               'with exact 95% confidence intervals, and carrier-frequency differences (group A minus group B) '
               'with Newcombe–Wilson 95% intervals. These intervals are pointwise, not multiplicity-adjusted. '
               'Apply Benjamini–Hochberg correction across all tested '
               'locus/haplotype hypotheses separately for the overall and pairwise families. These are exploratory '
               'unadjusted tests; covariate-adjusted or relatedness-aware effects require a separately validated design. '
               'Inference requires verified reference-call provenance, '
               f'at least {parameters["min_test_samples"]} resolved participants per group, '
               f'at least {parameters["min_resolved_fraction"]:.0%} resolved per group, '
               f'and at least {parameters["min_carriers"]} pooled carriers and noncarriers. '
               'A positive group-A-minus-group-B difference indicates more carriers in group A.')
    causal = ('No associated cohort haplotype is available for candidate-variant localization. '
              'The external GWAS lead table describes supplied external association evidence, not causal variants '
              'established in these three groups. The legacy SuSiE audit records why its exploratory rankings '
              'are not used as validated causal evidence. Resolving candidates requires reliable joint genotypes '
              'including reference calls, usable phase, harmonized full regional GWAS summary statistics with '
              'sample-size information, suitable ancestry-matched LD, and fine-mapping diagnostics. '
              'Functional evidence and independent replication would be needed to strengthen causal interpretation.')
    if candidates:
        causal = ('Candidate evidence is listed in candidate_variants.tsv with its provenance and limitations. '
                  'Statistical prioritization does not establish causality.')
    sections = [('Haplotype group comparisons', finding), ('Statistical methods', methods),
                ('Candidate causal variants', causal)]
    downloads = aggregate_files + evidence_files + ['comparison_status.tsv', 'comparison_results.json']
    extra = ''.join(f'<h2>{escape(title)}</h2><p>{escape(body)}</p>' for title, body in sections)
    extra += '<h2>Analysis tables</h2><ul>' + ''.join(
        f'<li><a href="{escape(name)}">{escape(name)}</a></li>' for name in downloads) + '</ul>'
    report = out / 'report.html'
    report.write_text(report.read_text().replace('<div class="provenance">', extra + '<div class="provenance">'))
    text_report = out / 'REPORT.txt'
    text_report.write_text(text_report.read_text() + '\n\n' + '\n\n'.join(
        title.upper() + '\n' + body for title, body in sections) + '\n')
    print(json.dumps(counts, indent=2))
    return counts


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', type=Path, required=True)
    parser.add_argument('--evidence-dir', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    publish(args.run_dir, args.evidence_dir, args.out)
