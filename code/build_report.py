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


def _figure(figure):
    # Image paths resolve from report.qmd; download links resolve from the HTML.
    png = os.path.relpath(figure['png'], PROJECT_ROOT)
    pdf = os.path.relpath(figure['pdf'], PROJECT_ROOT / 'results')
    title, caption = _text(figure['title']), _text(figure['caption'])
    guide = _figure_guide(figure)
    return (f'### {title}\n\n{guide}\n\n![{caption}](<{png}>){{width=100%}}\n\n'
            f'<a href="{html.escape(pdf, quote=True)}">Download PDF</a>\n')


def _report_intro(run_dir, synthetic, summary=None):
    purpose = ('This project compares patterns of nearby DNA variants between participant '
               'groups. A **haplotype** is one such pattern on one chromosome. Public GWAS '
               'signals help choose which DNA regions to examine.')
    if synthetic:
        status = ('**Real cohort analysis is pending.** The first figure uses real public '
                  'GWAS data. The later figures are a **synthetic example** made from '
                  'generated participants; they show how the analysis works.')
    elif run_dir:
        analysed = (summary or {}).get('locus_status_counts', {}).get('analysed', 0)
        status = ('The first figure uses public GWAS data. The cohort section reports '
                  'the completed analysis and its participant coverage.')
        if not analysed:
            status += (' **No region yielded an analysable haplotype panel.** '
                       'The cohort section explains why.')
    else:
        status = ('**Cohort analysis is pending.** This report shows the public GWAS data '
                  'and the status of the participant labels. No synthetic results are '
                  'substituted automatically.')
    return f'## What this report shows\n\n{purpose}\n\n{status}'


def _public_gwas_section(manifest):
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
    figures = [_figure(item) for item in manifest.get('figures', [])
               if item['scope'] == 'filtered_gwas']
    return '\n\n'.join([text, *figures])


def _participant_section(review, figures=()):
    heading = '## Participant labels and coverage'
    plots = [_figure(item) for item in figures if item['scope'] == 'participant_metadata']
    if not review:
        return '\n\n'.join([heading, 'No participant metadata review is available. '
                              'The report therefore cannot confirm the cohort group labels.', *plots])
    authoritative = review.get('authoritative_group_source', False)
    finalized = review.get('finalized_sample_table', False)
    included = review.get('included_sample_count')
    excluded = review.get('excluded_vcf_count')
    if authoritative:
        status = ('**Your supplied participant table is the authoritative group source.** '
                  'The category names are kept exactly as supplied.')
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
        parts.extend([explanation,
                      _table(['Group', 'People in supplied table'], sorted(supplied.items()))])
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


def _preflight_section(preflight):
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
    parts.append('This check does not calculate haplotype frequencies or compare participant groups. '
                 'It checks lead sites only, using the thresholds shown above.')
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


def _run_section(manifest, run_dir, synthetic, review=None):
    summary = manifest.get('run_summary')
    if not run_dir or not summary:
        return ('## Cohort figures are pending\n\n'
                'No completed cohort run was selected. Figures require a completed analysis '
                'using the participant assignments and coverage described above.')
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
    parts = [heading, explanation, _table(['Run group', 'Participants'], groups.items())]
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
    parts.extend(_figure(item) for item in manifest.get('figures', [])
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


def _technical_section(manifest, review, run_dir, config):
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
        parts.extend(['Existing ID-list counts, retained for source comparison:',
                      _table(['Group', 'People'], sorted(review['existing_category_counts'].items()))])
    command = ('.venv/bin/python code/build_report.py'
               f' --config {shlex.quote(str(resolve_project_path(config)))}'
               + (f' --run-dir {shlex.quote(str(run_dir))}' if run_dir else '') + ' --render')
    parts.append('To recreate this report, run the following from the project root. '
                 'The command records the exact selected run.\n\n```bash\n' + command + '\n```')
    parts.append('Editable source: `report.qmd`. Generated text: `results/report_content.qmd`. '
                 'Figure sources and captions: `figure_manifest.json`. PNGs are embedded in '
                 'the HTML; PDF downloads remain under `results/figures/`.')
    return '## Files and reproducibility\n\n' + _details('Show source files, settings and command', '\n\n'.join(parts))


def report_content(manifest, review_summary, run_dir, synthetic, config=DEFAULT_CONFIG, preflight_summary=None):
    """Assemble a plain-language report from counts, figures and saved reviews."""
    sections = [
        _report_intro(run_dir, synthetic, manifest.get('run_summary')),
        _public_gwas_section(manifest),
        _participant_section(review_summary, manifest.get('figures', [])),
        _preflight_section(preflight_summary),
        _run_section(manifest, run_dir, synthetic, review_summary),
        _workflow_section(),
        _technical_section(manifest, review_summary, run_dir, config),
    ]
    return '\n\n'.join(section for section in sections if section) + '\n'


def render_report():
    executable = shutil.which('quarto')
    local = PROJECT_ROOT / '.tools' / 'quarto' / 'bin' / 'quarto'
    if executable is None and local.is_file():
        executable = str(local)
    if executable is None:
        raise RuntimeError('Quarto is required for --render. Install Quarto or put it at '
                           '.tools/quarto/bin/quarto; the QMD content and figures are already prepared.')
    target = PROJECT_ROOT / 'results' / 'report.html'
    # Quarto can clean its output directory. Always render into an isolated
    # disposable directory, never into the directory holding analysis outputs.
    with tempfile.TemporaryDirectory(prefix='haplotype_report_render_') as staging:
        subprocess.run([executable, 'render', 'report.qmd', '--to', 'html',
                        '--output-dir', staging], cwd=PROJECT_ROOT, check=True)
        rendered = Path(staging) / 'report.html'
        if not rendered.is_file():
            raise RuntimeError(f'Quarto did not create the expected report: {rendered}')
        target.parent.mkdir(exist_ok=True)
        shutil.copy2(rendered, target)
    return target


def build_report(config=DEFAULT_CONFIG, run_dir=None, demo=False, render=False, max_loci=6):
    from make_figures import build_figures
    settings = load_config(config)
    selected, synthetic = select_run(settings, run_dir, demo)
    output = PROJECT_ROOT / 'results'
    output.mkdir(exist_ok=True)
    figure_parent = output / 'figures' / ('demo' if synthetic else 'cohort' if selected else 'gwas')
    figure_parent.mkdir(parents=True, exist_ok=True)
    figure_dir = Path(tempfile.mkdtemp(prefix='build_', dir=figure_parent))
    metadata_path = Path(settings['paths']['metadata'])
    manifest = build_figures(gwas_path=settings['paths']['gwas'], out_dir=figure_dir,
                             run_dir=selected, synthetic=synthetic,
                             p_threshold=settings['analysis']['p_threshold'], max_loci=max_loci,
                             metadata_path=metadata_path if metadata_path.is_file() else None)
    review_path = PROJECT_ROOT / 'data' / 'intermediate' / 'participant_review' / 'summary.json'
    review = json.loads(review_path.read_text()) if review_path.is_file() else None
    preflight_path = PROJECT_ROOT / 'data' / 'intermediate' / 'cohort_preflight' / 'summary.json'
    preflight = json.loads(preflight_path.read_text()) if preflight_path.is_file() else None
    content = output / 'report_content.qmd'
    content.write_text(report_content(manifest, review, selected, synthetic, config, preflight), encoding='utf-8')
    context = dict(created_at_utc=datetime.now(timezone.utc).isoformat(),
                   config=str(config), run_dir=str(selected) if selected else None,
                   synthetic=synthetic, figure_manifest=str(figure_dir / 'figure_manifest.json'),
                   participant_review=str(review_path) if review else None,
                   genotype_preflight=str(preflight_path) if preflight else None,
                   content=str(content), report_source=str(PROJECT_ROOT / 'report.qmd'))
    (output / 'report_context.json').write_text(json.dumps(context, indent=2) + '\n', encoding='utf-8')
    if render:
        context['html'] = str(render_report())
    return context


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', default=str(DEFAULT_CONFIG))
    choice = parser.add_mutually_exclusive_group()
    choice.add_argument('--run-dir', type=Path, help='Completed analysis directory; no fallback if missing')
    choice.add_argument('--demo', action='store_true', help='Use the latest completed synthetic demo explicitly')
    parser.add_argument('--max-loci', type=int, default=6, help='Maximum loci with detailed figures')
    parser.add_argument('--render', action='store_true', help='Render results/report.html with Quarto')
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
