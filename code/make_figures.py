#!/usr/bin/env python3
"""Export filtered-GWAS and explicitly selected run figures as PNG and PDF."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
import json
import math
from pathlib import Path
import re
import sys

from project_config import DEFAULT_CONFIG, PROJECT_ROOT, load_config


def neg_log10_p(value):
    """Transform a positive P-value without converting tiny P to binary float."""
    try:
        p = value if isinstance(value, Decimal) else Decimal(str(value))
    except (InvalidOperation, ValueError) as error:
        raise ValueError(f'Invalid P-value: {value!r}') from error
    if not p.is_finite() or not 0 < p <= 1:
        raise ValueError('P-values must be finite and within (0, 1]')
    result = float(-p.log10())
    if not math.isfinite(result):
        raise ValueError('Transformed P-value exceeds the plotting range')
    return result


def _count(value):
    """Read a table count without silently accepting fractions or negatives."""
    if isinstance(value, bool):
        raise ValueError('Counts must be nonnegative integers')
    result = int(value)
    if result < 0 or (isinstance(value, float) and value != result):
        raise ValueError('Counts must be nonnegative integers')
    return result


def ratio_or_na(numerator, denominator):
    """Return an observed fraction, preserving zero versus no denominator."""
    numerator, denominator = _count(numerator), _count(denominator)
    if numerator > denominator:
        raise ValueError('Numerator exceeds denominator')
    return numerator / denominator if denominator else math.nan


def prepare_phase_matrix(rows, loci, groups):
    """Return locus-by-group fractions of participants with resolved phase."""
    cells = {}
    for row in rows:
        key = (row['locus_id'], row['group'])
        if key in cells:
            raise ValueError(f'Duplicate phase-resolution row: {key}')
        resolved = _count(row['n_resolved'])
        total = _count(row['n_group'])
        cells[key] = (ratio_or_na(resolved, total), f'{resolved}/{total}')

    values, counts = [], []
    for locus in loci:
        row_values, row_counts = [], []
        for group in groups:
            fraction, label = cells.get((locus, group), (math.nan, 'NA'))
            row_values.append(fraction)
            row_counts.append(label)
        values.append(row_values)
        counts.append(row_counts)
    return {
        'loci': list(loci), 'groups': list(groups),
        'values': values, 'counts': counts,
    }


def prepare_haplotype_matrix(rows, locus_id, groups, max_haplotypes=24):
    """Return haplotype-by-group copy fractions among resolved chromosomes."""
    if type(max_haplotypes) is not int or max_haplotypes < 1:
        raise ValueError('max_haplotypes must be a positive integer')
    cells = {}
    pooled_copies = defaultdict(int)
    denominators = {group: None for group in groups}
    for row in rows:
        if row['locus_id'] != locus_id:
            continue
        haplotype, group = row['haplotype_id'], row['group']
        key = (haplotype, group)
        if key in cells:
            raise ValueError(f'Duplicate haplotype-frequency row: {key}')
        copies = _count(row['haplotype_copies'])
        chromosomes = _count(row['resolved_chromosomes'])
        cells[key] = ratio_or_na(copies, chromosomes)
        pooled_copies[haplotype] += copies
        if group in denominators:
            prior = denominators[group]
            if prior is not None and prior != chromosomes:
                raise ValueError(f'Inconsistent resolved-chromosome denominator: {locus_id}, {group}')
            denominators[group] = chromosomes
    ranked_haplotypes = sorted(pooled_copies, key=lambda haplotype: (-pooled_copies[haplotype], haplotype))
    haplotypes = ranked_haplotypes[:max_haplotypes]
    values = []
    for haplotype in haplotypes:
        values.append([cells.get((haplotype, group), math.nan) for group in groups])
    return {
        'haplotypes': haplotypes, 'groups': list(groups),
        'values': values,
        'total_haplotypes': len(ranked_haplotypes),
        'truncated': len(ranked_haplotypes) > max_haplotypes,
        'denominators': denominators,
    }


def _read_table(path, required):
    """Read one analysis table and check the columns used by its chart."""
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not reader.fieldnames or not set(required) <= set(reader.fieldnames):
            raise ValueError(f'{path}: missing required columns: {", ".join(required)}')
        return list(reader)


def _read_gwas(path, threshold):
    """Collect valid autosomal points, counting the inclusive cutoff separately."""
    points = []
    counts = {
        'gwas_rows': 0, 'plotted_rows': 0, 'at_or_below_threshold': 0,
        'above_threshold': 0, 'skipped_rows': 0,
    }
    with path.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        columns = reader.fieldnames or []
        chrom_column = next((name for name in ['Chromsome', 'Chromosome', 'CHROM', 'chrom']
                             if name in columns), None)
        if chrom_column is None or not {'Position', 'Pval'} <= set(columns):
            raise ValueError(f'{path}: GWAS needs chromosome, Position and Pval columns')
        for row in reader:
            counts['gwas_rows'] += 1
            try:
                chrom = int(re.sub(r'^chr', '', row[chrom_column], flags=re.I))
                position = int(row['Position'])
                p = Decimal(row['Pval'])
                if not 1 <= chrom <= 22 or position < 1:
                    raise ValueError('Outside positive autosomal coordinates')
                log_p = neg_log10_p(p)
            except (ValueError, InvalidOperation, TypeError):
                counts['skipped_rows'] += 1
                continue
            points.append((chrom, position, log_p))
            if p <= threshold:
                counts['at_or_below_threshold'] += 1
            else:
                counts['above_threshold'] += 1
    counts['plotted_rows'] = len(points)
    if not points:
        raise ValueError('No valid positive-P autosomal GWAS rows to plot')
    return points, counts


def _load_run(run_dir):
    """Load aggregate results only from the completed run explicitly requested."""
    if not (run_dir / 'COMPLETED.txt').is_file() or not (run_dir / 'run_summary.json').is_file():
        raise ValueError(f'Requested run is incomplete: {run_dir}; COMPLETED.txt and run_summary.json are required')
    summary = json.loads((run_dir / 'run_summary.json').read_text())
    if not isinstance(summary, dict) or not isinstance(summary.get('group_counts'), dict):
        raise ValueError('Run summary must contain group_counts')
    filenames = ['lead_loci.tsv', 'group_phase_qc.tsv', 'haplotype_group_frequencies.tsv']
    for name in filenames:
        if not (run_dir / name).is_file():
            raise ValueError(f'Requested run is missing {name}: {run_dir}')
    leads = _read_table(run_dir / filenames[0], ['locus_id', 'status'])
    phase = _read_table(run_dir / filenames[1], ['locus_id', 'group', 'n_group', 'n_resolved'])
    frequencies = _read_table(run_dir / filenames[2],
                              ['locus_id', 'haplotype_id', 'group', 'haplotype_copies', 'resolved_chromosomes'])
    return summary, leads, phase, frequencies


def _pyplot():
    """Import the plotting dependency only when a figure is being exported."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise RuntimeError('Figure export requires matplotlib. Install dependencies with '
                           'python3 -m pip install -r requirements.txt, or use .venv/bin/python.') from error
    return plt


def _style_axis(ax):
    for side in ['top', 'right']:
        ax.spines[side].set_visible(False)
    ax.spines['bottom'].set_color('#c3cbd5')
    ax.spines['left'].set_color('#c3cbd5')
    ax.tick_params(colors='#475569', labelsize=9, length=3)
    ax.set_axisbelow(True)


def _save(plt, fig, out, ident, title, caption, scope):
    """Export both formats and return their manifest entry without overwriting."""
    png, pdf = out / f'{ident}.png', out / f'{ident}.pdf'
    if png.exists() or pdf.exists():
        plt.close(fig)
        raise ValueError(f'Refusing to overwrite figure files for {ident}')
    fig.savefig(png, dpi=180, facecolor='white')
    fig.savefig(pdf, facecolor='white', metadata={'Title': title, 'Subject': caption})
    plt.close(fig)
    return {'id': ident, 'title': title, 'caption': caption,
            'png': str(png), 'pdf': str(pdf), 'scope': scope}


def _overview(plt, points, counts, threshold, out):
    """Show all supplied GWAS points, with a second panel for the lower range."""
    by_chrom = defaultdict(list)
    for chrom, position, log_p in points:
        by_chrom[chrom].append((position, log_p))
    chroms = sorted(by_chrom)
    chromosome_ends = {chrom: max(position for position, _ in by_chrom[chrom]) for chrom in chroms}
    gap = sum(chromosome_ends.values()) * 0.004
    ticks, series = [], []
    offset = 0
    for index, chrom in enumerate(chroms):
        values = by_chrom[chrom]
        positions, logs = zip(*values)
        x = [position + offset for position in positions]
        ticks.append(offset + (min(positions) + max(positions)) / 2)
        series.append((x, logs, ['#285a88', '#219d98'][index % 2]))
        offset += chromosome_ends[chrom] + gap
    maximum = max(point[2] for point in points)
    cutoff = neg_log10_p(threshold) if threshold else None
    detail_limit = min(50, max(12, maximum * 1.06))
    fig, axes = plt.subplots(2, 1, figsize=(13, 8), sharex=True,
                             gridspec_kw={'height_ratios': [1, 1.05], 'hspace': 0.18})
    fig.subplots_adjust(left=0.085, right=0.975, top=0.84, bottom=0.16)
    fig.text(0.085, 0.95, 'GWAS variants by chromosome',
             fontsize=19, weight='bold', color='#172f46')
    fig.text(0.085, 0.907, f"{counts['plotted_rows']:,} supplied variants  |  "
             f"{counts['at_or_below_threshold']:,} with P ≤ {threshold}  |  "
             'Filtered input only', fontsize=10, color='#526579')
    for ax in axes:
        for x, logs, color in series:
            ax.scatter(x, logs, s=4, alpha=0.6, c=color, linewidths=0, rasterized=True)
        if cutoff is not None:
            ax.axhline(cutoff, color='#b45332', linewidth=1.1, linestyle=(0, (4, 3)),
                       label=f'P ≤ {threshold}')
        ax.set_ylabel('−log₁₀(P)\n(higher = smaller P)', fontsize=10, color='#243b53')
        ax.grid(axis='y', color='#e8edf2', linewidth=0.7)
        ax.set_xlim(0, offset)
        _style_axis(ax)
    axes[0].set_ylim(0, max(12, maximum * 1.12))
    axes[0].text(0.01, 0.90, 'All values in the supplied file', transform=axes[0].transAxes,
                 fontsize=10, color='#40566d')
    axes[0].legend(loc='upper right', frameon=False, fontsize=9)
    axes[1].set_ylim(0, detail_limit)
    axes[1].text(0.01, 0.90, f'Detail: 0–{detail_limit:g}; higher values shown above',
                 transform=axes[1].transAxes, fontsize=10, color='#40566d',
                 bbox={'facecolor': 'white', 'alpha': 0.85, 'edgecolor': 'none'})
    axes[1].set_xticks(ticks, [str(chrom) for chrom in chroms])
    axes[1].set_xlabel('Chromosome · each dot is one variant at its position in the supplied file', fontsize=10, labelpad=10)
    reading_note = 'Reading the vertical scale: 8 means P = 10⁻⁸. Each increase of 1 means a tenfold smaller P-value.'
    data_note = (f'Only filtered rows are shown. Maximum −log₁₀(P): {maximum:.2f}. '
                 f'Invalid rows or rows outside chromosomes 1–22 omitted: {counts["skipped_rows"]:,}.')
    fig.text(0.085, 0.065, reading_note + '\n' + data_note,
             fontsize=9, color='#63758a', linespacing=1.6, va='top')
    caption = (f'Each dot is a variant from the supplied filtered genome-wide association study (GWAS): '
               f'{counts["plotted_rows"]:,} valid autosomal rows, '
               f'{counts["at_or_below_threshold"]:,} at the inclusive P ≤ {threshold} cutoff. '
               'Alternating colors separate chromosomes. Higher points have smaller P-values: '
               '8 on the vertical scale means P = 10^-8, and each increase of 1 means a tenfold smaller P-value. '
               'The upper panel shows the full range; the lower panel zooms in. '
               'Decimal logarithms preserve extremely small P-values. Only filtered rows are shown, '
               'so this figure cannot assess the P-value distribution of the full GWAS.')
    return _save(plt, fig, out, 'filtered_gwas_overview', 'Filtered GWAS overview', caption, 'filtered_gwas')


def _percentage_label(value):
    """Keep small nonzero values visibly different from an observed zero."""
    if math.isnan(value):
        return 'NA'
    if value == 0:
        return '0%'
    if value == 1:
        return '100%'
    percent = 100 * value
    if percent < 0.1:
        return '<0.1%'
    if percent > 99.9:
        return '>99.9%'
    return f'{percent:.1f}'.rstrip('0').rstrip('.') + '%'


def _matrix_figure(plt, matrix, row_labels, column_labels, *, title, subtitle, footer,
                   counts=None, synthetic=False):
    """Draw a percentage table, with a separate banner for simulated data."""
    height = max(5.5, min(18, 3.5 + 0.47 * len(row_labels)))
    fig, ax = plt.subplots(figsize=(max(11.5, 3 + len(column_labels) * 1.4), height))
    title_y = 0.875 if synthetic else 0.945
    subtitle_y = title_y - 0.065
    fig.subplots_adjust(left=0.19, right=0.89, top=subtitle_y - 0.10, bottom=0.30)
    if synthetic:
        fig.text(0.065, 0.965, 'SYNTHETIC DEMO · simulated participants',
                 fontsize=11, weight='bold', color='#9a3412', va='top',
                 bbox={'boxstyle': 'round,pad=0.35', 'facecolor': '#fff1e6', 'edgecolor': '#c45508'})
    fig.text(0.065, title_y, title, fontsize=17, weight='bold', color='#172f46')
    fig.text(0.065, subtitle_y, subtitle, fontsize=10, color='#526579', va='top', linespacing=1.4)
    cmap = plt.get_cmap('YlGnBu').with_extremes(bad='#e5e7eb')
    image = ax.imshow(matrix, cmap=cmap, vmin=0, vmax=1, aspect='auto')
    ax.set_xticks(range(len(column_labels)), column_labels, rotation=15, ha='right')
    ax.set_yticks(range(len(row_labels)), row_labels)
    ax.tick_params(length=0, labelsize=10, pad=7)
    for spine in ax.spines.values():
        spine.set_visible(False)
    for row, values in enumerate(matrix):
        for column, value in enumerate(values):
            label = _percentage_label(value)
            if counts is not None and not math.isnan(value):
                label = f'{counts[row][column]}  ·  {label}'
            ax.text(column, row, label, ha='center', va='center', fontsize=9,
                    color='white' if not math.isnan(value) and value > 0.55 else '#253749')
    colorbar = fig.colorbar(image, ax=ax, fraction=0.045, pad=0.04, ticks=[0, .25, .5, .75, 1])
    colorbar.ax.set_yticklabels(['0%', '25%', '50%', '75%', '100%'])
    colorbar.outline.set_visible(False)
    colorbar.ax.tick_params(labelsize=9, length=0)
    fig.text(0.065, 0.11, footer, fontsize=9, color='#63758a', linespacing=1.5, va='top')
    return fig


def _phase_figure(plt, data, total_loci, out, synthetic):
    """Explain how many participants have two usable chromosome patterns."""
    title = 'Participants with resolved haplotypes'
    selection = f'Showing {len(data["loci"])} of {total_loci} analyzed regions in lead-table order.'
    footer = (
        'Resolved = the two patterns of nearby variants could be determined for a participant.\n'
        'Cells show resolved participants / all participants in that group. Gray NA means no estimate.\n'
        + selection
    )
    fig = _matrix_figure(
        plt, data['values'], data['loci'], data['groups'], title=title,
        subtitle='How often the nearby-variant patterns on both chromosomes could be determined.',
        footer=footer, counts=data['counts'], synthetic=synthetic,
    )
    prefix = 'SYNTHETIC DEMO · ' if synthetic else ''
    caption = prefix + (
        'A haplotype is a pattern of nearby variants on one chromosome. '
        'Each cell shows how many participants had both chromosome patterns resolved, '
        'divided by all participants in that group (n_resolved / n_group). '
        'Participants with unresolved patterns remain in the denominator. '
        'Darker cells mean a higher resolved share; gray NA means no estimate. ' + selection
    )
    scope = 'synthetic_demo' if synthetic else 'cohort_run'
    return _save(plt, fig, out, 'phase_resolution', prefix + title, caption, scope)


def _haplotype_figure(plt, locus, data, out, synthetic):
    """Explain the copy frequency of each observed nearby-variant pattern."""
    title = f'{locus}: haplotype frequencies'
    total = data['total_haplotypes']
    if data['truncated']:
        shown = len(data['haplotypes'])
        selection = (f'Top {shown} of {total} patterns by total observed copies across groups. '
                     'Percentages keep the original denominators.')
    else:
        selection = f'All {total} observed patterns; sorted by total observed copies across groups.'

    denominators = []
    for group in data['groups']:
        denominator = data['denominators'][group]
        denominators.append(f'{group}: {denominator if denominator is not None else "NA"}')
    subtitle = ('Resolved chromosome copies (two per resolved participant):\n'
                + '  ·  '.join(denominators))
    footer = (
        'Haplotype = a pattern of nearby variants on one chromosome. Percentages describe chromosome copies.\n'
        '0% means the pattern was absent among resolved copies. Gray NA means no estimate.\n'
        + selection
    )
    fig = _matrix_figure(
        plt, data['values'], data['haplotypes'], data['groups'],
        title=title, subtitle=subtitle, footer=footer, synthetic=synthetic,
    )
    prefix = 'SYNTHETIC DEMO · ' if synthetic else ''
    caption = prefix + (
        f'{locus}: each row is a haplotype, a pattern of nearby variants on one chromosome. '
        'Each cell is copies of that pattern divided by all resolved chromosome copies in the group '
        '(haplotype_copies / resolved_chromosomes). Each resolved participant contributes two copies. '
        'These percentages describe chromosome copies rather than the share of people carrying the pattern. '
        'Darker cells mean a more frequent pattern; 0% is an observed zero, while gray NA means no estimate. '
        + selection
    )
    safe_locus = re.sub(r'[^A-Za-z0-9_-]', '_', locus)
    scope = 'synthetic_demo' if synthetic else 'cohort_run'
    return _save(plt, fig, out, f'haplotype_frequencies_{safe_locus}', prefix + title, caption, scope)


def _category_figure(plt, category_counts, out):
    """Show the supplied table's category counts, without implying VCF matching."""
    preferred = ['T2D', 'Prediabetes', 'NORM']
    groups = [group for group in preferred if group in category_counts]
    groups.extend(sorted(set(category_counts) - set(preferred)))
    counts = [category_counts[group] for group in groups]
    total = sum(counts)
    fig, ax = plt.subplots(figsize=(10, 5.5))
    fig.subplots_adjust(left=0.10, right=0.96, top=0.77, bottom=0.24)
    fig.text(0.10, 0.93, 'Participant categories from your table',
             fontsize=19, weight='bold', color='#172f46')
    fig.text(0.10, 0.85, f'{total} participants · exact supplied Category labels',
             fontsize=11, color='#526579')
    bars = ax.bar(groups, counts, color=['#b76d38', '#219d98', '#285a88'], width=0.55)
    ax.bar_label(bars, labels=[str(count) for count in counts], padding=5,
                 fontsize=13, weight='bold')
    ax.set_ylim(0, max(counts) * 1.22)
    ax.set_ylabel('Participants in the supplied table', fontsize=10)
    ax.yaxis.get_major_locator().set_params(integer=True)
    ax.grid(axis='y', color='#e8edf2', linewidth=0.7)
    _style_axis(ax)
    ax.tick_params(axis='x', labelsize=12, pad=8)
    fig.text(0.10, 0.12,
             'These counts describe the category source, before matching names to the VCF.\n'
             'They do not measure genotype completeness or haplotype frequencies.',
             fontsize=10, color='#63758a', linespacing=1.5, va='top')
    caption = ('Counts of unique participant IDs in the supplied Category table: '
               + ', '.join(f'{group}: {category_counts[group]}' for group in groups)
               + '. Category labels are preserved. These are source-table counts; '
               'they do not imply that every ID has been matched to the VCF or analysed.')
    return _save(plt, fig, out, 'participant_category_counts',
                 'Participant categories from your table', caption, 'participant_metadata')


def build_figures(*, gwas_path, out_dir, run_dir=None, synthetic=False, p_threshold='1e-8', max_loci=6,
                  metadata_path=None):
    """Plot the GWAS, an optional participant table, and an optional completed run."""
    if type(max_loci) is not int or max_loci < 1:
        raise ValueError('max_loci must be a positive integer')
    try:
        threshold = Decimal(str(p_threshold))
    except InvalidOperation as error:
        raise ValueError('p_threshold must be a decimal within [0,1]') from error
    if not threshold.is_finite() or not 0 <= threshold <= 1:
        raise ValueError('p_threshold must be within [0,1]')
    gwas_path = Path(gwas_path).resolve()
    out = Path(out_dir).resolve()
    selected_run = Path(run_dir).resolve() if run_dir is not None else None
    run_data = _load_run(selected_run) if selected_run is not None else None
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise ValueError(f'Figure output must be a new or empty directory: {out}')
    points, counts = _read_gwas(gwas_path, threshold)
    category_counts = None
    if metadata_path is not None:
        from make_sample_template import read_metadata_groups
        metadata_path = Path(metadata_path).resolve()
        category_counts = dict(Counter(read_metadata_groups(metadata_path).values()))
        if not category_counts:
            raise ValueError('The participant category table is empty')
    if run_data is None:
        summary, leads, phase, frequencies = None, [], [], []
    else:
        summary, leads, phase, frequencies = run_data
    analysed = [row['locus_id'] for row in leads if row['status'] == 'analysed']
    groups = list(summary['group_counts']) if summary is not None else []
    # Keep exported tables readable without changing any frequency denominator.
    phase_loci = analysed[:40]
    selected_loci = analysed[:max_loci]
    phase_data = prepare_phase_matrix(phase, phase_loci, groups) if phase_loci else None
    haplotype_data = [(locus, prepare_haplotype_matrix(frequencies, locus, groups)) for locus in selected_loci]
    if analysed and not groups:
        raise ValueError('Analyzed loci require group counts in the run summary')
    for locus, data in haplotype_data:
        if not data['haplotypes']:
            raise ValueError(f'No haplotype-frequency rows for analyzed locus {locus}')
    plt = _pyplot()
    out.mkdir(parents=True, exist_ok=True)
    figures = []
    with plt.rc_context({'font.family': 'DejaVu Sans', 'pdf.fonttype': 42, 'axes.titleweight': 'bold'}):
        figures.append(_overview(plt, points, counts, threshold, out))
        if category_counts is not None:
            figures.append(_category_figure(plt, category_counts, out))
        if phase_data is not None:
            figures.append(_phase_figure(plt, phase_data, len(analysed), out, synthetic))
        for locus, data in haplotype_data:
            figures.append(_haplotype_figure(plt, locus, data, out, synthetic))

    counts.update({
        'analysed_loci': len(analysed), 'phase_loci_shown': len(phase_loci),
        'haplotype_loci_shown': len(selected_loci), 'figures': len(figures),
    })
    run_files = {}
    if selected_run is not None:
        for name in ['COMPLETED.txt', 'run_summary.json', 'lead_loci.tsv',
                     'group_phase_qc.tsv', 'haplotype_group_frequencies.tsv']:
            run_files[name] = str(selected_run / name)
    source_paths = [gwas_path] + [Path(path) for path in run_files.values()]
    if metadata_path is not None:
        source_paths.append(metadata_path)
        counts['participant_category_counts'] = category_counts
        counts['participant_source_rows'] = sum(category_counts.values())
    source_stats = {}
    for path in source_paths:
        stats = path.stat()
        source_stats[str(path)] = {'bytes': stats.st_size, 'mtime_ns': stats.st_mtime_ns}
    manifest = {
        'schema_version': 1, 'created_at_utc': datetime.now(timezone.utc).isoformat(),
        'figures': figures,
        'sources': {
            'gwas': str(gwas_path), 'run_dir': str(selected_run) if selected_run else None,
            'run_files': run_files,
            'metadata': str(metadata_path) if metadata_path is not None else None,
        },
        'source_stats': source_stats,
        'synthetic': bool(synthetic), 'p_threshold': str(threshold), 'max_loci': max_loci,
        'counts': counts, 'run_summary': summary,
    }
    (out / 'figure_manifest.json').write_text(json.dumps(manifest, indent=2, allow_nan=False) + '\n')
    return manifest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', default=str(DEFAULT_CONFIG))
    parser.add_argument('--gwas', help='Filtered GWAS TSV; default from project settings')
    parser.add_argument('--metadata', help='Optional participant Category TSV for source-count figures')
    parser.add_argument('--run-dir', help='Explicit completed analysis directory; omitted means GWAS only')
    parser.add_argument('--synthetic', action='store_true', help='Label every run figure SYNTHETIC DEMO')
    parser.add_argument('--out', default=str(PROJECT_ROOT / 'results' / 'figures'), help='New or empty figure directory')
    parser.add_argument('--p-threshold', help='Inclusive cutoff; default from project settings')
    parser.add_argument('--max-loci', type=int, default=6, help='Maximum analyzed loci with haplotype heatmaps')
    args = parser.parse_args(argv)
    try:
        settings = load_config(args.config)
        manifest = build_figures(gwas_path=args.gwas or settings['paths']['gwas'], out_dir=args.out,
                                 run_dir=args.run_dir, synthetic=args.synthetic,
                                 p_threshold=args.p_threshold or settings['analysis']['p_threshold'],
                                 max_loci=args.max_loci, metadata_path=args.metadata)
    except (OSError, ValueError, RuntimeError, csv.Error) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        return 2
    print(f'Exported {len(manifest["figures"])} figures as PNG and PDF to {Path(args.out).resolve()}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
