#!/usr/bin/env python3
"""Prepare local figures and Quarto report content; optionally render the HTML."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import html
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile

from project_config import DEFAULT_CONFIG, PROJECT_ROOT, load_config, resolve_project_path
from group_labels import display_group, ordered_groups


RUN_TABLES = (
    'COMPLETED.txt', 'lead_loci.tsv', 'group_phase_qc.tsv', 'haplotype_group_frequencies.tsv',
    'selected_markers.tsv', 'sample_haplotypes.tsv',
)


def _completed_summary(folder):
    """Require a completion record and the tables consumed by the report."""
    folder = Path(folder)
    try:
        summary = json.loads((folder / 'run_summary.json').read_text())
        completed = datetime.fromisoformat(summary['completed_at_utc'].replace('Z', '+00:00'))
        if completed.tzinfo is None:
            raise ValueError('completion time needs a timezone')
        if not isinstance(summary['group_counts'], dict) or type(summary['selected_loci']) is not int:
            raise ValueError('invalid group counts or locus count')
        missing = [name for name in RUN_TABLES if not (folder / name).is_file()]
        if missing:
            raise ValueError(f'missing completion marker or tables: {", ".join(missing)}')
    except (OSError, ValueError, KeyError, TypeError, AttributeError) as error:
        raise ValueError(f'No completed analysis at {folder}: {error}') from error
    return summary


def _inside(path, folder):
    try:
        Path(path).resolve().relative_to(Path(folder).resolve())
        return True
    except ValueError:
        return False


def select_run(settings, run_dir=None, demo=False):
    """Choose an explicit/latest demo or configured cohort run, without fallback."""
    if run_dir is not None and demo:
        raise ValueError('Choose --run-dir or --demo, not both.')
    paths = settings['paths']
    if demo:
        completed = []
        parent = Path(paths['demo_results'])
        for candidate in parent.iterdir() if parent.is_dir() else []:
            if not candidate.is_dir():
                continue
            try:
                summary = _completed_summary(candidate)
            except ValueError:
                continue
            stamp = datetime.fromisoformat(summary['completed_at_utc'].replace('Z', '+00:00'))
            completed.append((stamp, str(candidate), candidate))
        if not completed:
            raise ValueError('No completed synthetic demo found. Run bash code/run_demo.sh first.')
        return max(completed)[2].resolve(), True
    selected = Path(run_dir if run_dir is not None else paths['out']).expanduser().resolve()
    try:
        summary = _completed_summary(selected)
    except ValueError:
        if run_dir is not None:
            raise
        return None, False
    source_vcf = summary.get('parameters', {}).get('vcf')
    synthetic = _inside(selected, paths['demo_results']) or bool(
        source_vcf and _inside(source_vcf, paths['demo_data']))
    return selected, synthetic


def _text(value):
    """Escape source labels for Markdown tables and raw HTML contexts."""
    return html.escape(str(value)).replace('|', '&#124;').replace('\n', ' ')


def _table(headers, rows):
    return '\n'.join([
        '| ' + ' | '.join(headers) + ' |',
        '| ' + ' | '.join('---' for _ in headers) + ' |',
        *('| ' + ' | '.join(_text(cell) for cell in row) + ' |' for row in rows),
    ])


def _count_text(value):
    """Make counts easy to read without changing missing values to zero."""
    return f'{value:,}' if isinstance(value, int) else _text(value)


def _details(title, content):
    return f'<details>\n<summary>{_text(title)}</summary>\n\n{content}\n\n</details>'


def _figure_guide(figure):
    """Explain the reading task for each plot type."""
    ident = figure['id']
    if ident == 'cohort_comparison_overview':
        return ('Read the workflow and participant counts first. The lead-site panel '
                'shows the share of participants with usable genotypes at each selected GWAS lead. '
                'Pooled values combine all included participants; colored points show each group. '
                'The reference line marks the run\'s genotype-completeness requirement. '
                '**These percentages measure available genotype data, not haplotype frequencies '
                'or diabetes risk.**\n\n'
                'The arrangement is inspired by the study overview in '
                '[Figueroa et al. (2024), Figure 1]'
                '(https://www.nature.com/articles/s41588-024-01719-5/figures/1). '
                'All plotted cohort values come from this analysis; no published patient '
                'patterns are reproduced.')
    if ident == 'regional_marker_comparison':
        return ('Each row is a participant group. The left column shows usable genotypes '
                'divided by all participants in that group; the right column shows ALT allele '
                'copies divided by observed chromosome copies at each SNP. The shared horizontal '
                'scale gives physical positions around the lead, marked by a vertical line. '
                'The caption records how this region was selected. Read allele frequency '
                'alongside genotype completeness: frequency is conditional on observed genotypes, '
                'and missing calls are excluded. No estimate is drawn when its denominator is zero. '
                '**These are single-variant summaries, not phased chromosome patterns or '
                'haplotype frequencies.** Unconnected points describe separate SNPs; matching '
                'values across groups do not establish shared ancestry.\n\n'
                'The common genomic scale is inspired by the local-variant comparison in '
                '[Figueroa et al. (2024), Figure 1d]'
                '(https://www.nature.com/articles/s41588-024-01719-5/figures/1). '
                'This figure uses our cohort\'s observed calls and does not reproduce the '
                'paper\'s expansion-bearing chromosome patterns.')
    if figure.get('scope') == 'lead_genotype_qc':
        return ('Each row is a GWAS lead variant and each column is a group. Cells show '
                'the share of matched participants with a genotype meeting the depth and quality '
                'thresholds. The denominator is the matched participant count in that group. '
                'These are genotype-completeness estimates, not haplotype frequencies.')
    if ident == 'participant_category_counts':
        return ('Each bar counts people in one category of the supplied participant table. '
                'These are source-table counts; the people included in a completed genetic '
                'analysis are reported separately.')
    if ident == 'filtered_gwas_overview':
        return ('Each dot is a variant, arranged by chromosome and position. Higher dots '
                'have smaller P-values and stronger evidence of association in the public '
                'study. The dashed line marks the chosen cutoff. The top panel shows the '
                'full range; the bottom panel zooms in so smaller signals are easier to see.')
    if ident == 'phase_resolution':
        return ('Rows are DNA regions and columns are groups. A cell such as **3/4 (75%)** '
                'means that usable haplotype pairs were identified for three of four people. '
                'Lower percentages mean less complete information for that comparison.')
    if ident.startswith('haplotype_frequencies_'):
        return ('Rows are haplotype patterns and columns are groups. Read across a row to '
                'compare the same pattern. **25%** means one quarter of the resolved '
                'chromosome copies have that pattern; each resolved person contributes '
                'two copies. **0%** means none were observed; **NA** means no estimate is available.')
    return ''


def _figure(figure, report_dir=None):
    # Image paths resolve from report.qmd; download links resolve from the HTML.
    png = os.path.relpath(figure['png'], report_dir or PROJECT_ROOT)
    pdf = os.path.relpath(figure['pdf'], report_dir or PROJECT_ROOT / 'results')
    title, caption = _text(figure['title']), _text(figure['caption'])
    guide = _figure_guide(figure)
    return (f'### {title}\n\n{guide}\n\n![{caption}](<{png}>){{width=100%}}\n\n'
            f'<a href="{html.escape(pdf, quote=True)}">Download PDF</a>\n')


def _report_intro(run_dir, synthetic, summary=None, review=None, preflight=None):
    purpose = ('This project compares patterns of nearby DNA variants between participant '
               'groups. A **haplotype** is one such pattern on one chromosome. Public GWAS '
               'signals help choose which DNA regions to examine.')
    if synthetic:
        status = ('**Real cohort haplotype analysis is pending.** The first figure uses real public '
                  'GWAS data. Only the haplotype and phase-resolution panels marked '
                  '**SYNTHETIC DEMO** use generated participants. Participant-category '
                  'and genotype-completeness figures, when present, use their separately '
                  'recorded sources. The demo panels illustrate how the analysis works.')
    elif run_dir:
        analysed = (summary or {}).get('locus_status_counts', {}).get('analysed', 0)
        status = ('The first figure uses public GWAS data. The cohort section reports '
                  'the completed analysis and its participant coverage.')
        if not analysed:
            status += (' **No region yielded an analysable haplotype panel.** '
                       'The cohort section explains why.')
    else:
        status = '**Cohort haplotype analysis is pending.** This report shows the public GWAS data.'
        if (review or {}).get('authoritative_group_source'):
            status += (' Participant categories are confirmed; the supplied group counts '
                       'and VCF identifier coverage are shown below.')
            pending = review.get('pending_vcf_id_matches', review.get('unmatched_metadata_ids'))
            if pending:
                status += f' **{_count_text(pending)} identifier matches still await confirmation.**'
        else:
            status += ' The participant-label review is summarized below.'
        if preflight:
            status += (' Real genotype-completeness checks describe the available calls; '
                       'they do not calculate haplotype frequencies.')
        status += ' **This report contains no synthetic demo panels.**'
    return f'## What this report shows\n\n{purpose}\n\n{status}'


def _public_gwas_section(manifest, report_dir=None):
    counts = manifest.get('counts', {})
    total = _count_text(counts.get('gwas_rows', 'not reported'))
    passing = _count_text(counts.get('at_or_below_threshold', 'not reported'))
    threshold = _text(manifest.get('p_threshold', '1e-8'))
    filename = _text(Path(manifest['sources']['gwas']).name)
    text = (f'## Public GWAS: variants to examine\n\n'
            f'`{filename}` contains **{total} variant rows**. Of these, **{passing}** '
            f'meet the selected threshold, **P ≤ {threshold}**. A smaller P-value means '
            'stronger evidence of association in the GWAS; it does not show that a variant '
            'causes diabetes.\n\n'
            'The file was already filtered before this analysis. This figure shows those '
            'selected variants, so it is not a complete view of all variants across the genome.')
    figures = [_figure(item, report_dir) for item in manifest.get('figures', [])
               if item['scope'] == 'filtered_gwas']
    return '\n\n'.join([text, *figures])


def _participant_section(review, figures=(), report_dir=None):
    heading = '## Participant labels and coverage'
    plots = [_figure(item, report_dir) for item in figures if item['scope'] == 'participant_metadata']
    if not review:
        return '\n\n'.join([heading, 'No participant metadata review is available. '
                              'The report therefore cannot confirm the cohort group labels.', *plots])
    authoritative = review.get('authoritative_group_source', False)
    finalized = review.get('finalized_sample_table', False)
    included = review.get('included_sample_count')
    excluded = review.get('excluded_vcf_count')
    if authoritative:
        status = ('**Your supplied participant table is the authoritative group source.** '
                  'Source category values are preserved in the inputs; their report labels '
                  'are shown alongside them below.')
        group_column = review.get('metadata_group_column')
        if group_column:
            status += f' Group assignments come from the `{_text(group_column)}` column.'
        changed = review.get('previous_category_change_count')
        if changed is not None:
            status += (f' Compared with the previous participant table, '
                       f'**{_count_text(changed)} participant group assignments changed**.')
        label_only = review.get('category_label_only_change_count')
        if label_only:
            status += (f' A further **{_count_text(label_only)} category labels changed in name only**, '
                       'without changing the displayed phenotype group.')
        conflicts = review.get('comparison_counts', {}).get('category_conflict', 0)
        if conflicts:
            status += (f' The **{_count_text(conflicts)} historical group-label conflicts** '
                       'with the older ID lists are resolved in favour of this table.')
        if finalized:
            status += ' The sample table is finalized for the included participants.'
        else:
            status += (' Matching all table identifiers to VCF sample identifiers remains '
                       'incomplete; the group-source choice is settled.')
            pending = review.get('pending_vcf_id_matches', review.get('unmatched_metadata_ids'))
            if pending is not None:
                status += f' **{_count_text(pending)} VCF identifier matches await confirmation.**'
    elif finalized:
        status = 'The saved metadata review records a finalized sample table.'
    else:
        conflicts = _count_text(review.get('comparison_counts', {}).get('category_conflict', 'not reported'))
        unmatched = _count_text(review.get('unmatched_metadata_ids', 'not reported'))
        status = (f'**Final cohort labels remain unresolved.** The metadata review found '
                  f'**{conflicts} group-label conflicts** between the supplied sources and '
                  f'**{unmatched} identifiers** that do not exactly match the VCF. '
                  'Confirm the correct labels and identifiers before running the cohort comparison.')
    if included is not None:
        status += f' **{_count_text(included)} participants are included** in the prepared sample table.'
    if excluded:
        status += (f' **{_count_text(excluded)} VCF samples are excluded** from that table. '
                   'They must be matched explicitly before they can be included.')
    supplied = review.get('authoritative_category_counts') or review.get('pasted_category_counts', {})
    parts = [heading, status]
    if supplied:
        explanation = ('These counts describe the complete supplied participant table. '
                       'Counts for any selected completed run appear separately below.'
                       if authoritative or finalized else
                       'The supplied participant table lists the following groups. These '
                       'counts describe that source and are not verified analysis assignments.')
        exact_counts = review.get('exact_match_category_counts', {})
        alias_counts = review.get('verified_alias_category_counts')
        matched_counts = review.get('matched_category_counts', {})
        show_aliases = isinstance(alias_counts, dict)
        coverage = []
        for group in ordered_groups(supplied):
            total, exact = supplied[group], exact_counts.get(group)
            alias = alias_counts.get(group, 0) if show_aliases else 0
            matched = matched_counts.get(group)
            if isinstance(exact, int) and isinstance(alias, int):
                computed = exact + alias
                if matched is not None and show_aliases and matched != computed:
                    raise ValueError(f'Participant coverage disagrees with exact and verified alias counts: {group}')
                if matched is None:
                    matched = computed
            if isinstance(total, int) and isinstance(matched, int) and not 0 <= matched <= total:
                raise ValueError(f'Participant matched count exceeds source category count: {group}')
            if show_aliases and (not isinstance(alias, int) or alias < 0):
                raise ValueError(f'Invalid verified alias count: {group}')
            pending = total - matched if isinstance(total, int) and isinstance(matched, int) else 'not reported'
            row = [group, display_group(group), total, exact if exact is not None else 'not reported']
            if show_aliases:
                row.append(alias)
            coverage.append([*row, pending])
        headers = ['Source category', 'Display group', 'People in table', 'Exact VCF matches']
        if show_aliases:
            headers.append('Verified alias matches')
            explanation += (' Verified alias matches use an explicitly confirmed identifier map; '
                            'they are separate from exact identifier matches.')
        headers.append('Pending ID matches')
        numeric_columns = range(2, len(headers))
        if coverage and all(isinstance(row[index], int) for row in coverage for index in numeric_columns):
            coverage.append(['Total', '', *[sum(row[index] for row in coverage) for index in numeric_columns]])
        parts.extend([explanation,
                      _table(headers, coverage)])
    return '\n\n'.join([*parts, *plots])


def _locus_status_table(statuses):
    explanations = {
        'negative_liftover_allele_orientation_unconfirmed': 'Lead allele orientation after negative-strand liftover is unconfirmed',
        'palindromic_lead_orientation_unconfirmed': 'Lead allele orientation is unconfirmed for an A/T or C/G lead',
        'no_unique_PASS_sequence_allele_match_at_lead': 'No eligible VCF record matched the lead position and alleles',
        'multiple_matching_lead_records': 'More than one VCF record matched the lead',
        'lead_below_min_call_rate': 'Too few participants had a usable lead genotype at the required quality',
        'fewer_than_two_eligible_markers': 'Fewer than two markers passed the panel criteria',
    }
    return _table(['Why the region was skipped', 'Regions'], [
        [explanations.get(reason, reason.replace('_', ' ')), count]
        for reason, count in statuses.items() if reason != 'analysed'
    ])


def _preflight_section(preflight, figures=(), report_dir=None):
    """Describe lead-site coverage without attributing it to phenotype groups."""
    if not preflight:
        return ''
    parameters = preflight.get('parameters', {})
    samples = _count_text(preflight.get('n_samples', 'not reported'))
    leads = _count_text(preflight.get('n_leads', 'not reported'))
    passing = _count_text(preflight.get('leads_passing_call_rate', 'not reported'))
    maximum = _count_text(preflight.get('max_high_quality_called', 'not reported'))
    threshold = parameters.get('min_call_rate')
    threshold_text = f'{100 * threshold:.1f}%' if isinstance(threshold, (int, float)) else 'not reported'
    maximum_rate = preflight.get('max_high_quality_call_rate')
    maximum_text = f'{100 * maximum_rate:.1f}%' if isinstance(maximum_rate, (int, float)) else 'not reported'
    parts = [
        '## Lead-genotype completeness check',
        f'A separate indexed check examined **{leads} lead-variant sites** across '
        f'**all {samples} VCF samples**, without using participant group labels. '
        'This checks whether enough usable genotypes are present before attempting haplotype comparisons.',
        _table(['Check', 'Result'], [
            ['Required genotype depth (DP)', parameters.get('min_dp', 'not reported')],
            ['Required genotype quality (GQ)', parameters.get('min_gq', 'not reported')],
            ['Required share of usable genotypes at a lead', threshold_text],
            ['Leads meeting that share', f'{passing} of {leads}'],
            ['Highest usable-genotype count at any checked lead', f'{maximum} of {samples} ({maximum_text})'],
            ['Recorded reference-homozygote genotypes at these sites', preflight.get('raw_reference_genotypes', 'not reported')],
        ]),
    ]
    if preflight.get('leads_passing_call_rate') == 0:
        parts.append('**No checked lead met the required genotype-completeness threshold.** '
                     f'With all {samples} samples included and these thresholds applied, '
                     'these leads cannot form an eligible haplotype panel.')
    if preflight.get('raw_reference_genotypes') == 0:
        parts.append('No reference-homozygote calls (`0/0` or `0|0`) were recorded at these '
                     'lead sites. Missing genotypes cannot be treated as reference calls; '
                     'the VCF generation and reference-call completeness need review.')
    parts.append('This check does not calculate haplotype frequencies or test phenotype associations. '
                 'It checks lead sites only, using the thresholds shown above.')
    group_qc = preflight.get('group_qc')
    if group_qc:
        matched_groups = group_qc.get('matched_category_counts', {})
        matched_count = _count_text(group_qc.get('matched_sample_count', 'not reported'))
        unassigned = _count_text(group_qc.get('unmatched_vcf_count', 'not reported'))
        has_aliases = group_qc.get('scope') == 'exact_and_verified_id_matches'
        if has_aliases:
            exact = _count_text(group_qc.get('exact_match_sample_count', 'not reported'))
            aliases = _count_text(group_qc.get('verified_alias_sample_count', 'not reported'))
            heading = '### Group completeness among exact and verified ID matches'
            explanation = (f'The group-specific check uses **{matched_count} participants**: '
                           f'**{exact} exact metadata–VCF ID matches** and '
                           f'**{aliases} verified alias matches** from the confirmed identifier map. '
                           f'**{unassigned} VCF samples remain unassigned.** '
                           'This check measures genotype completeness and does not produce haplotype results.')
            denominator = ('Each group percentage uses its combined exact and verified alias count '
                           'as the denominator. The all-sample check above also includes any unassigned samples.')
        else:
            heading = '### Group completeness among exact ID matches'
            explanation = (f'The group-specific check uses **{matched_count} participants** whose supplied '
                           f'metadata IDs exactly match the VCF. **{unassigned} VCF samples remain unassigned.** '
                           'This QC-only grouping does not finalize an analysis subset or produce haplotype results.')
            denominator = ('Each group percentage uses its exact-match count as the denominator. '
                           'The all-sample check above also includes the unassigned samples.')
        parts.extend([
            heading, explanation,
            _table(['Display group', 'Source category',
                    'Matches used for QC' if has_aliases else 'Exact matches used for QC'], [
                [display_group(group), group, matched_groups[group]]
                for group in ordered_groups(matched_groups)
            ]),
            denominator,
        ])
        parts.extend(_figure(item, report_dir) for item in figures if item['scope'] == 'lead_genotype_qc')
    if preflight.get('inputs_stable_during_read') is False:
        parts.append('The source files changed during the check. Repeat it before using these counts.')
    inputs = preflight.get('inputs', {})
    details = [['Check prepared (UTC)', preflight.get('created_at_utc', 'not reported')]]
    details.extend([name, item.get('path', 'not reported')]
                   for name, item in inputs.items() if isinstance(item, dict))
    parts.append(_details('Genotype-check sources', _table(['Source', 'Value'], details)))
    return '\n\n'.join(parts)


def _uses_authoritative_sample_table(summary, review):
    """Verify saved input content, rather than infer provenance from current labels."""
    if not review or not review.get('authoritative_group_source') or not review.get('finalized_sample_table'):
        return False
    expected = review.get('finalized_sample_table_sha256')
    if not isinstance(expected, str) or len(expected) != 64:
        return False
    if any(character not in '0123456789abcdefABCDEF' for character in expected):
        return False
    sample_path = summary.get('parameters', {}).get('samples')
    recorded = summary.get('small_input_sha256', {}).get(sample_path)
    return isinstance(recorded, str) and recorded.lower() == expected.lower()


def _run_section(manifest, run_dir, synthetic, review=None, report_dir=None):
    summary = manifest.get('run_summary')
    if not run_dir or not summary:
        return ('## Cohort haplotype results are unavailable\n\n'
                'No completed cohort haplotype run was selected. Genotype-completeness checks '
                'are shown separately above when available. Haplotype frequency estimates '
                'require a completed analysis using the participant assignments and coverage '
                'described above.')
    groups = summary.get('group_counts', {})
    participants = _count_text(sum(groups.values()))
    statuses = summary.get('locus_status_counts', {})
    analysed = statuses.get('analysed', 0)
    loci = _count_text(analysed)
    if synthetic:
        heading = '## Synthetic example: how to read the analysis'
        explanation = (f'This example contains **{participants} generated participants** and '
                       f'**{loci} analysed DNA regions**. Every figure in this section uses '
                       '**synthetic data**, not results from your participants.')
    else:
        heading = '## Selected cohort analysis'
        explanation = (f'This completed run contains **{participants} participants** and '
                       f'**{loci} analysed DNA regions**.')
        if (review or {}).get('authoritative_group_source'):
            if _uses_authoritative_sample_table(summary, review):
                explanation += (' The recorded sample-table checksum matches the finalized '
                                'authoritative table, so these group assignments use that table.')
            else:
                explanation += (' Your supplied table is authoritative, but this selected run '
                                'has not been verified against its finalized sample table.')
        elif not (review or {}).get('finalized_sample_table'):
            explanation += ' See the participant-label review before interpreting these groups.'
    parts = [heading, explanation, _table(['Display group', 'Source group', 'Participants'], [
        [display_group(group), group, groups[group]] for group in ordered_groups(groups)
    ])]
    if not analysed:
        parts.append('**No haplotype frequency estimates were produced for this run.** '
                     'None of the selected regions met all requirements for a haplotype panel. '
                     'This outcome does not show that the groups have identical haplotypes.')
        if statuses:
            parts.append(_locus_status_table(statuses))
        parts.append('Review `lead_loci.tsv` and `run_summary.json` for the recorded criteria and reasons.')
    else:
        parts.append('**Read completeness first, then compare patterns.** Missing or unresolved '
                     'haplotypes are excluded from the frequency calculation; they are not treated '
                     'as absent patterns. Different amounts of usable data can affect the comparison.')
    figure_ids = {item['id'] for item in manifest.get('figures', [])}
    if figure_ids & {'cohort_comparison_overview', 'regional_marker_comparison'}:
        parts.append('The comparison figures below use the participant groups and regional '
                     'variant tables saved by this completed run. They can show observed genotypes '
                     'and their completeness even when no haplotype panel was resolved. '
                     'Preparing these figures does not load or rerun a separate VCF precheck.')
    if ('regional_marker_comparison' in figure_ids and
            summary.get('scan_counts', {}).get('high_quality_homozygous_reference_genotypes') == 0):
        parts.append('**The run recorded no high-quality homozygous-reference (`0/0`) calls '
                     'in the selected windows.** Observed allele frequencies can be biased '
                     'by incomplete reference-call representation; these plots do not test '
                     'group associations.')
    parts.extend(_figure(item, report_dir) for item in manifest.get('figures', [])
                 if item['scope'] in {'cohort_run', 'synthetic_demo'})
    if not summary.get('haplotype_carrier_tests_performed', 0):
        parts.append('No haplotype association tests were run.')
    return '\n\n'.join(parts)


def _workflow_section():
    return ('## How the code works\n\n'
            '1. Keep GWAS variants that meet the chosen P-value cutoff.\n'
            '2. Group nearby signals and choose a lead variant for each region.\n'
            '3. Read nearby genotypes from the VCF and identify haplotype pairs where '
            'the available data allow it.\n'
            '4. Compare pattern frequencies across groups, alongside the amount of usable data.\n\n'
            'Both genetic inputs must use the same verified genome build and compatible '
            'allele representations. The workflow compares local patterns; it does not '
            'perform statistical fine-mapping or establish causal variants.\n\n'
            '### Short glossary\n\n' + _table(['Term', 'Meaning'], [
                ['GWAS', 'A study that tests whether genetic variants are associated with a trait.'],
                ['Variant / allele', 'A DNA difference / one version at that position.'],
                ['Haplotype', 'A pattern of nearby variants on one chromosome.'],
                ['Locus', 'A DNA region examined by the analysis.'],
                ['Resolved', 'The genotype and phase information allow the two haplotype copies to be identified.'],
                ['VCF', 'The file containing variant records and participant genotypes.'],
            ]))


def _technical_section(manifest, review, run_dir, config, preflight_path=None,
                       no_preflight=False, output_dir=None):
    """Keep exact provenance and settings available without crowding the reading."""
    summary = manifest.get('run_summary') or {}
    counts = manifest.get('counts', {})
    rows = [
        ['Report prepared (UTC)', manifest.get('created_at_utc', 'not reported')],
        ['GWAS source', manifest['sources']['gwas']],
        ['Analysis directory', run_dir or 'No completed run selected'],
        ['Settings', resolve_project_path(config)],
        ['Rows plotted / skipped', f'{counts.get("plotted_rows", "not reported")} / {counts.get("skipped_rows", "not reported")}'],
        ['Rows above the P-value cutoff', counts.get('above_threshold', 'not reported')],
    ]
    if summary:
        rows.extend([
            ['Run completed (UTC)', summary.get('completed_at_utc', 'not reported')],
            ['Declared genome build', summary.get('genome_build_user_declared', 'not reported')],
            ['Selected regions', summary.get('selected_loci', 'not reported')],
            ['Haplotype carrier tests', summary.get('haplotype_carrier_tests_performed', 0)],
        ])
    parts = [_table(['Source detail', 'Value'], rows)]
    if summary:
        parameters = summary.get('parameters', {})
        parts.append(_table(['Recorded setting', 'Value'], [
            [key, parameters.get(key, 'not reported')] for key in
            ('p_threshold', 'lead_gap', 'top_loci', 'flank', 'max_markers',
             'min_dp', 'min_gq', 'min_call_rate')]))
        parts.append(_table(['Region status', 'Count'], summary.get('locus_status_counts', {}).items()))
    if review and review.get('existing_category_counts'):
        counts_by_group = review['existing_category_counts']
        parts.extend(['Existing ID-list counts, retained for source comparison:',
                      _table(['Display group', 'Source category', 'People'], [
                          [display_group(group), group, counts_by_group[group]]
                          for group in ordered_groups(counts_by_group)
                      ])])
    command = ('.venv/bin/python code/build_report.py'
               f' --config {shlex.quote(str(resolve_project_path(config)))}'
               + (f' --run-dir {shlex.quote(str(run_dir))}' if run_dir else '')
               + (f' --preflight {shlex.quote(str(preflight_path))}' if preflight_path else '')
               + (' --no-preflight' if no_preflight else '')
               + (f' --output-dir {shlex.quote(str(output_dir))}' if output_dir is not None else '')
               + ' --render')
    parts.append('To recreate this report, run the following from the project root. '
                 'The command records the exact selected run.\n\n```bash\n' + command + '\n```')
    location = '' if output_dir is not None else 'results/'
    parts.append(f'Editable source: `report.qmd`. Generated text: `{location}report_content.qmd`. '
                 'Figure sources and captions: `figure_manifest.json`. PNGs are embedded in '
                 f'the HTML; PDF downloads remain under `{location}figures/`.')
    if no_preflight:
        parts.append('No separate genotype-preflight summary was loaded for this report. '
                     'The completed analysis retains its recorded genotype filters.')
    return '## Files and reproducibility\n\n' + _details('Show source files, settings and command', '\n\n'.join(parts))


def report_content(manifest, review_summary, run_dir, synthetic, config=DEFAULT_CONFIG,
                   preflight_summary=None, preflight_path=None, no_preflight=False, output_dir=None):
    """Assemble a plain-language report from counts, figures and saved reviews."""
    sections = [
        _report_intro(run_dir, synthetic, manifest.get('run_summary'), review_summary, preflight_summary),
        _public_gwas_section(manifest, output_dir),
        _participant_section(review_summary, manifest.get('figures', []), output_dir),
        _preflight_section(preflight_summary, manifest.get('figures', []), output_dir),
        _run_section(manifest, run_dir, synthetic, review_summary, output_dir),
        _workflow_section(),
        _technical_section(manifest, review_summary, run_dir, config, preflight_path,
                           no_preflight, output_dir),
    ]
    return '\n\n'.join(section for section in sections if section) + '\n'


def render_report(report_source=None, output_dir=None):
    executable = shutil.which('quarto')
    local = PROJECT_ROOT / '.tools' / 'quarto' / 'bin' / 'quarto'
    if executable is None and local.is_file():
        executable = str(local)
    if executable is None:
        raise RuntimeError('Quarto is required for --render. Install Quarto or put it at '
                           '.tools/quarto/bin/quarto; the QMD content and figures are already prepared.')
    source = Path(report_source) if report_source is not None else PROJECT_ROOT / 'report.qmd'
    output = Path(output_dir) if output_dir is not None else PROJECT_ROOT / 'results'
    target = output / 'report.html'
    # Quarto can clean its output directory. Always render into an isolated
    # disposable directory, never into the directory holding analysis outputs.
    with tempfile.TemporaryDirectory(prefix='haplotype_report_render_') as staging:
        subprocess.run([executable, 'render', source.name, '--to', 'html',
                        '--output-dir', staging], cwd=source.parent, check=True)
        rendered = Path(staging) / 'report.html'
        if not rendered.is_file():
            raise RuntimeError(f'Quarto did not create the expected report: {rendered}')
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(rendered, target)
    return target


def build_report(config=DEFAULT_CONFIG, run_dir=None, demo=False, render=False, max_loci=6,
                 preflight=None, no_preflight=False, output_dir=None):
    from make_figures import build_figures
    if no_preflight and preflight is not None:
        raise ValueError('Choose --preflight or --no-preflight, not both.')
    settings = load_config(config)
    selected, synthetic = select_run(settings, run_dir, demo)
    preflight_path = None if no_preflight else (
        Path(preflight).expanduser().resolve() if preflight is not None else
        PROJECT_ROOT / 'data' / 'intermediate' / 'cohort_preflight' / 'summary.json')
    if preflight is not None and not preflight_path.is_file():
        raise ValueError(f'Genotype preflight summary does not exist: {preflight_path}')
    preflight_summary = (json.loads(preflight_path.read_text())
                         if preflight_path is not None and preflight_path.is_file() else None)
    output = Path(output_dir).expanduser().resolve() if output_dir is not None else PROJECT_ROOT / 'results'
    source = PROJECT_ROOT / 'report.qmd'
    custom_source = None
    if output_dir is not None:
        include = '{{< include results/report_content.qmd >}}'
        custom_source = source.read_text(encoding='utf-8')
        if custom_source.count(include) != 1:
            raise ValueError('Report template must contain exactly one results/report_content.qmd include.')
        custom_source = custom_source.replace(include, '{{< include report_content.qmd >}}')
        source = output / 'report.qmd'
        if source.resolve() == (PROJECT_ROOT / 'report.qmd').resolve():
            raise ValueError('A custom output directory must not overwrite the project report template.')
    output.mkdir(parents=True, exist_ok=True)
    figure_parent = output / 'figures' / ('demo' if synthetic else 'cohort' if selected else 'gwas')
    figure_parent.mkdir(parents=True, exist_ok=True)
    figure_dir = Path(tempfile.mkdtemp(prefix='build_', dir=figure_parent))
    metadata_path = Path(settings['paths']['metadata'])
    manifest = build_figures(gwas_path=settings['paths']['gwas'], out_dir=figure_dir,
                             run_dir=selected, synthetic=synthetic,
                             p_threshold=settings['analysis']['p_threshold'], max_loci=max_loci,
                             metadata_path=metadata_path if metadata_path.is_file() else None,
                             preflight_path=preflight_path if preflight_summary and preflight_summary.get('group_qc') else None)
    review_path = PROJECT_ROOT / 'data' / 'intermediate' / 'participant_review' / 'summary.json'
    review = json.loads(review_path.read_text()) if review_path.is_file() else None
    content = output / 'report_content.qmd'
    content.write_text(report_content(manifest, review, selected, synthetic, config,
                                      preflight_summary, preflight_path if preflight_summary else None,
                                      no_preflight, output if output_dir is not None else None), encoding='utf-8')
    if custom_source is not None:
        source.write_text(custom_source, encoding='utf-8')
    context = dict(created_at_utc=datetime.now(timezone.utc).isoformat(),
                   config=str(config), run_dir=str(selected) if selected else None,
                   synthetic=synthetic, figure_manifest=str(figure_dir / 'figure_manifest.json'),
                   participant_review=str(review_path) if review else None,
                   genotype_preflight=str(preflight_path) if preflight_summary else None,
                   preflight_disabled=bool(no_preflight),
                   content=str(content), report_source=str(source))
    (output / 'report_context.json').write_text(json.dumps(context, indent=2) + '\n', encoding='utf-8')
    if render:
        context['html'] = str(render_report(source, output))
    return context


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', default=str(DEFAULT_CONFIG))
    choice = parser.add_mutually_exclusive_group()
    choice.add_argument('--run-dir', type=Path, help='Completed analysis directory; no fallback if missing')
    choice.add_argument('--demo', action='store_true', help='Use the latest completed synthetic demo explicitly')
    parser.add_argument('--max-loci', type=int, default=6, help='Maximum loci with detailed figures')
    preflight_choice = parser.add_mutually_exclusive_group()
    preflight_choice.add_argument('--preflight', type=Path,
                                 help='Use this genotype-completeness summary JSON instead of the saved default')
    preflight_choice.add_argument('--no-preflight', action='store_true',
                                 help='Do not load any saved genotype-preflight summary or group-QC figure')
    parser.add_argument('--output-dir', type=Path,
                        help='Save the report, content, context and figures in this directory (default: results)')
    parser.add_argument('--render', action='store_true', help='Render report.html in the report output directory with Quarto')
    args = parser.parse_args(argv)
    if args.max_loci < 1:
        parser.error('--max-loci must be positive')
    try:
        result = build_report(**vars(args))
    except (ValueError, OSError, RuntimeError, subprocess.CalledProcessError) as error:
        parser.error(str(error))
    print(json.dumps(result, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
