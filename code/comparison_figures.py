"""Paper-inspired layouts drawn exclusively from completed-run aggregate tables."""
from __future__ import annotations

from collections import defaultdict
import csv
import gzip
import math

from group_labels import display_group, ordered_groups


PAPER_FIGURE = 'https://www.nature.com/articles/s41588-024-01719-5/figures/1'
COLORS = {'T2D': '#c66a32', 'Prediabetes': '#218d88', 'Control': '#3569a1', 'NORM': '#3569a1'}


def prepare_comparison_data(rows, summary, leads):
    """Keep completeness and conditional allele-copy fractions distinct.

    Only count-derived values are used. Validate every allele denominator,
    including multiallelic records; use biallelic SNPs for regional displays.
    Lead choice is independent of phenotype contrasts: maximum pooled coverage,
    with ties resolved by the existing GWAS lead-table order.
    """
    groups = ordered_groups(summary['group_counts'])
    totals = {g: int(summary['group_counts'][g]) for g in groups}
    variants = {}
    lead_ids = {lead['locus_id'] for lead in leads}
    for row in rows:
        ident, group = row['locus_id'], row['group']
        if ident not in lead_ids or group not in totals:
            raise ValueError('Regional table has a locus or group absent from the completed run')
        key = ident, row['variant_id']
        n, called, chromosomes, copies, index = (
            int(row[k]) for k in ('n_group', 'n_called', 'called_chromosomes', 'allele_copies', 'allele_index'))
        if n != totals[group] or not 0 <= called <= n or chromosomes != 2 * called or not 0 <= copies <= chromosomes:
            raise ValueError('Regional allele counts or genotype denominators are inconsistent')
        definition = (row['chrom'], int(row['pos']), row['ref'], row['alt'], row['type'],
                      int(row['distance_to_lead_bp']), row['is_exact_lead_match'] == '1')
        item = variants.setdefault(key, {'definition': definition, 'groups': {}})
        if definition != item['definition']:
            raise ValueError('Inconsistent regional variant definition')
        cell = item['groups'].setdefault(group, {'n_called': called, 'chromosomes': chromosomes, 'copies': {}})
        if cell['n_called'] != called or cell['chromosomes'] != chromosomes or index in cell['copies']:
            raise ValueError('Duplicate alleles or inconsistent denominators in regional table')
        cell['copies'][index] = copies

    by_locus = defaultdict(list)
    exact = defaultdict(list)
    for (ident, variant_id), item in variants.items():
        chrom, pos, ref, alt, kind, distance, is_lead = item['definition']
        alleles = alt.split(',')
        if set(item['groups']) != set(groups):
            raise ValueError('Regional variant is missing a participant group')
        for cell in item['groups'].values():
            if set(cell['copies']) != set(range(len(alleles) + 1)) or sum(cell['copies'].values()) != cell['chromosomes']:
                raise ValueError('Regional allele copies do not sum to called chromosomes')
        item.update({'variant_id': variant_id, 'chrom': chrom, 'pos': pos, 'ref': ref, 'alt': alt,
                     'distance_kb': distance / 1000, 'is_lead': is_lead,
                     'biallelic_snp': kind == 'SNP' and len(ref) == 1 and len(alt) == 1})
        for group, cell in item['groups'].items():
            cell['completeness'] = cell['n_called'] / totals[group] if totals[group] else math.nan
            cell['alt_fraction'] = (cell['copies'][1] / cell['chromosomes']
                                    if item['biallelic_snp'] and cell['chromosomes'] else math.nan)
        by_locus[ident].append(item)
        if is_lead:
            exact[ident].append(item)

    total = sum(totals.values())
    coverage, candidates = [], []
    for index, lead in enumerate(leads):
        matches = exact[lead['locus_id']]
        unique = matches[0] if len(matches) == 1 else None
        n_called = sum(c['n_called'] for c in unique['groups'].values()) if unique else None
        row = {'lead': lead, 'variant': unique, 'n_called': n_called,
               'pooled': n_called / total if unique and total else math.nan,
               'fractions': [unique['groups'][g]['completeness'] if unique else math.nan for g in groups]}
        coverage.append(row)
        if unique and unique['biallelic_snp'] and n_called and total:
            candidates.append((row['pooled'], -index, row))
    example = max(candidates, key=lambda x: x[:2])[2] if candidates else None
    regional = sorted(by_locus[example['lead']['locus_id']], key=lambda item: (item['pos'], item['variant_id'])) if example else []
    return {'groups': groups, 'totals': totals, 'total': total, 'coverage': coverage,
            'example': example, 'regional': [r for r in regional if r['biallelic_snp']],
            'regional_total_records': len(regional),
            'regional_excluded_records': sum(not r['biallelic_snp'] for r in regional)}


def load_comparison_data(run_dir, summary, leads):
    path = run_dir / 'regional_variant_frequencies.tsv.gz'
    if not path.is_file():
        return None
    with gzip.open(path, 'rt', encoding='utf-8', newline='') as handle:
        return prepare_comparison_data(csv.DictReader(handle, delimiter='\t'), summary, leads)


def write_plot_tables(data, out):
    """Export exact plotted numerators and denominators alongside the figures."""
    coverage_path = out / 'cohort_comparison_plot_data.tsv'
    regional_path = out / 'regional_marker_plot_data.tsv'
    fields = ['locus_id', 'chrom', 'lead_pos', 'group', 'n_called', 'n_group', 'completeness']
    with coverage_path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t')
        writer.writeheader()
        for row in data['coverage']:
            for group, value in zip(data['groups'], row['fractions']):
                called = row['variant']['groups'][group]['n_called'] if row['variant'] else 'NA'
                writer.writerow({**{k: row['lead'][k] for k in ['locus_id', 'chrom', 'lead_pos']},
                                 'group': group, 'n_called': called, 'n_group': data['totals'][group],
                                 'completeness': value if math.isfinite(value) else 'NA'})
    fields = ['locus_id', 'variant_id', 'pos', 'ref', 'alt', 'distance_kb', 'is_lead', 'group',
              'n_called', 'n_group', 'completeness', 'alt_copies', 'called_chromosomes', 'observed_alt_fraction']
    with regional_path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t')
        writer.writeheader()
        for variant in data['regional']:
            for group, cell in variant['groups'].items():
                writer.writerow({'locus_id': data['example']['lead']['locus_id'],
                    **{k: variant[k] for k in ['variant_id', 'pos', 'ref', 'alt', 'distance_kb', 'is_lead']},
                    'group': group, 'n_called': cell['n_called'], 'n_group': data['totals'][group],
                    'completeness': cell['completeness'] if math.isfinite(cell['completeness']) else 'NA',
                    'alt_copies': cell['copies'][1], 'called_chromosomes': cell['chromosomes'],
                    'observed_alt_fraction': cell['alt_fraction'] if math.isfinite(cell['alt_fraction']) else 'NA'})
    return {'cohort_comparison_overview': str(coverage_path), 'regional_marker_comparison': str(regional_path)}


def _panel(ax, letter, title, pad=14):
    ax.set_title(f'{letter}   {title}', loc='left', fontsize=12, pad=pad)


def overview_figure(plt, data, summary, out, save, style, synthetic=False):
    """Lettered workflow, group counts and lead completeness; no inferred phase."""
    from matplotlib.lines import Line2D
    from matplotlib.patches import FancyBboxPatch

    fig = plt.figure(figsize=(14, 11.5))
    grid = fig.add_gridspec(2, 2, height_ratios=[1, 5], width_ratios=[1, 2.1],
                           left=.07, right=.965, bottom=.14, top=.87, hspace=.35, wspace=.52)
    flow, group_ax, coverage_ax = fig.add_subplot(grid[0, :]), fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])
    title = 'Cohort comparison: available genotype data'
    prefix = 'SYNTHETIC DEMO · ' if synthetic else ''
    fig.text(.055, .956, prefix + title, fontsize=20, weight='bold', color='#172f46')
    fig.text(.055, .922, 'Completed analysis · identical group definitions in every panel', fontsize=11, color='#526579')
    flow.set_axis_off()
    _panel(flow, 'a', 'Analysis workflow and recorded outcome')
    params = summary['parameters']
    analysed = summary.get('locus_status_counts', {}).get('analysed', 0)
    labels = [f"Inputs\n{data['total']} participants + public GWAS",
              f"Regions\n{len(data['coverage'])} GWAS-selected windows",
              f"Internal requirements\nDP ≥ {params['min_dp']:g}; GQ ≥ {params['min_gq']:g}; calls ≥ {params['min_call_rate']:.0%}",
              f'Haplotype panels\n{analysed} available']
    for i, (label, color) in enumerate(zip(labels, ['#e8f0f7', '#e5f2ef', '#fbefdf', '#edf0f4'])):
        x = .008 + i * .249
        flow.add_patch(FancyBboxPatch((x, .13), .229, .68, boxstyle='round,pad=0.006,rounding_size=0.02',
                                     facecolor=color, edgecolor='#b7c2cd', transform=flow.transAxes))
        flow.text(x + .115, .47, label, transform=flow.transAxes, ha='center', va='center', fontsize=10, linespacing=1.7)
        if i < 3:
            flow.annotate('', xy=(x+.246, .47), xytext=(x+.233, .47), xycoords='axes fraction',
                          arrowprops={'arrowstyle': '->', 'color': '#526579'})
    _panel(group_ax, 'b', 'Participants in this run')
    groups = data['groups']
    bars = group_ax.bar(range(len(groups)), [data['totals'][g] for g in groups],
                        color=[COLORS.get(g, '#63758a') for g in groups], width=.65)
    group_ax.bar_label(bars, padding=5, fontsize=13, weight='bold')
    group_ax.set_xticks(range(len(groups)), [display_group(g) for g in groups], rotation=15, ha='right')
    group_ax.set_ylabel('Participants')
    group_ax.set_ylim(0, max(data['totals'].values()) * 1.3)
    group_ax.yaxis.get_major_locator().set_params(integer=True)
    group_ax.grid(axis='y', color='#e8edf2')
    style(group_ax)
    # Reserve the lower part of this column for a direct explanation of outcomes.
    pos = group_ax.get_position()
    group_ax.set_position([pos.x0, pos.y0 + pos.height * .47, pos.width, pos.height * .53])
    states = summary.get('locus_status_counts', {})
    reason_labels = {'lead_below_min_call_rate': 'Incomplete lead genotypes',
                     'palindromic_lead_orientation_unconfirmed': 'Allele orientation unconfirmed',
                     'analysed': 'Usable haplotype panels'}
    lines = ['Recorded region outcomes'] + [f'{n}  {reason_labels.get(k, k.replace("_", " "))}' for k, n in states.items()]
    lines += ['', f'{summary.get("haplotype_carrier_tests_performed", 0)} haplotype association tests']
    if not analysed:
        lines += ['', 'Group differences have not', 'been estimated.']
    fig.text(pos.x0, pos.y0 + pos.height * .34, '\n'.join(lines), fontsize=10, linespacing=1.8, va='top', color='#40566d')

    _panel(coverage_ax, 'c', 'Usable genotypes at each GWAS lead', pad=36)
    offsets = [(i - (len(groups)-1)/2) * .2 for i in range(len(groups))]
    for i, g in enumerate(groups):
        coverage_ax.scatter([100*r['fractions'][i] for r in data['coverage']],
                            [j+offsets[i] for j in range(len(data['coverage']))], s=21,
                            color=COLORS.get(g, '#63758a'), zorder=3)
    coverage_ax.scatter([100*r['pooled'] for r in data['coverage']], range(len(data['coverage'])),
                        color='#273849', marker='D', s=16, zorder=4)
    cutoff = params['min_call_rate'] * 100
    coverage_ax.axvline(cutoff, color='#a53d33', ls='--', lw=1)
    labels = [f"{r['lead']['locus_id']} · {r['lead']['chrom'].removeprefix('chr')}:{int(r['lead']['lead_pos']):,}" for r in data['coverage']]
    coverage_ax.set_yticks(range(len(labels)), labels, fontsize=8)
    for i, row in enumerate(data['coverage']):
        if row['n_called'] is None:
            coverage_ax.text(2, i, 'NA: no unique exact lead', va='center', fontsize=8, color='#63758a')
    coverage_ax.set_ylim(len(labels)-.4, -.6)
    coverage_ax.set_xlim(-2, 103)
    coverage_ax.set_xticks([0, 25, 50, 75, 100])
    coverage_ax.set_xlabel('Participants with usable lead genotypes (%)', labelpad=10)
    coverage_ax.grid(axis='x', color='#e8edf2')
    style(coverage_ax)
    handles = [Line2D([], [], marker='o', ls='', color=COLORS.get(g, '#63758a'), label=display_group(g), markersize=5) for g in groups]
    handles += [Line2D([], [], marker='D', ls='', color='#273849', label='All participants', markersize=4)]
    coverage_ax.legend(handles=handles, loc='lower left', bbox_to_anchor=(-.01, 1.005), frameon=False,
                       ncol=4, fontsize=8, handletextpad=.3, columnspacing=.8)
    fig.text(.055, .069, f'Dashed line: {cutoff:g}% pooled completeness requirement. Region labels include chromosome and {summary["genome_build_user_declared"]} position.\n'
             'Coverage is calculated from the saved analysis tables. It measures available genotypes, not haplotype sharing or diabetes risk.\n'
             'Panel layout inspired by Figueroa et al., Nature Genetics (2024), Fig. 1; all plotted cohort values come from this run.',
             fontsize=9, color='#526579', linespacing=1.6, va='top')
    caption = (prefix + 'a, Inputs, region selection and current outcome. b, Participants actually included in the completed run. '
               'c, Usable lead genotypes divided by participants in each group; black diamonds use the pooled denominator. '
               f'The dashed {cutoff:g}% line is the pooled call-rate requirement, not a test of group differences. '
               'Only unique exact lead matches are plotted; otherwise coverage is NA. '
               'These are completeness measurements, not haplotype frequencies. '
               f'Layout inspiration: {PAPER_FIGURE}. Values are derived from the selected run, without loading preflight results.')
    return save(plt, fig, out, 'cohort_comparison_overview', prefix + title, caption,
                'synthetic_demo' if synthetic else 'cohort_run')


def regional_figure(plt, data, summary, out, save, style, synthetic=False):
    """Matched physical-position panels; unconnected points never imply phase."""
    if data['example'] is None:
        return None
    example, variants, groups = data['example'], data['regional'], data['groups']
    lead, lead_variant = example['lead'], example['variant']
    fig, axes = plt.subplots(len(groups), 2, figsize=(14, 3.2 + len(groups)*2.05), squeeze=False,
                             sharex=True, sharey=True)
    fig.subplots_adjust(left=.11, right=.96, top=.78, bottom=.19, wspace=.19, hspace=.45)
    title = f"{lead['locus_id']}: observed variants around the GWAS lead"
    prefix = 'SYNTHETIC DEMO · ' if synthetic else ''
    fig.text(.055, .957, prefix + title, fontsize=20, weight='bold', color='#172f46')
    subtitle = (f"{lead['chrom']}:{int(lead['lead_pos']):,}  {lead_variant['ref']}>{lead_variant['alt']} · "
                f"best pooled lead coverage: {example['n_called']}/{data['total']} ({example['pooled']:.1%})\n"
                f"{len(variants)} biallelic SNPs shown; {data['regional_excluded_records']} other records omitted. "
                'Points are individual sites; phase is not represented.')
    fig.text(.055, .91, subtitle, fontsize=10, color='#526579', linespacing=1.7, va='top')
    start_kb = (int(lead['start'])-int(lead['lead_pos']))/1000
    end_kb = (int(lead['end'])-int(lead['lead_pos']))/1000
    column_titles = ['a   Genotype completeness', 'b   Observed ALT-copy fraction']
    for j, column_title in enumerate(column_titles):
        axes[0, j].set_title(column_title, loc='left', fontsize=13, pad=27)
    for i, group in enumerate(groups):
        color = COLORS.get(group, '#63758a')
        cells = [v['groups'][group] for v in variants]
        lead_cell = lead_variant['groups'][group]
        for j, key in enumerate(['completeness', 'alt_fraction']):
            ax = axes[i, j]
            usable = [(v['distance_kb'], cell[key]) for v, cell in zip(variants, cells) if math.isfinite(cell[key])]
            ax.scatter([x for x, y in usable], [100*y for x, y in usable], s=11, color=color, alpha=.55,
                       linewidths=0, rasterized=True)
            # The exact lead is drawn last, larger and outlined, and numerators are printed.
            if math.isfinite(lead_cell[key]):
                ax.scatter([0], [100*lead_cell[key]], s=64, marker='D', color=color, edgecolor='#172f46', linewidth=.9, zorder=6)
            ax.axvline(0, color='#a53d33', lw=.9, ls='--', zorder=0)
            ax.set_xlim(start_kb, end_kb)
            ax.set_ylim(-4, 104)
            ax.set_yticks([0, 25, 50, 75, 100])
            ax.grid(axis='y', color='#e8edf2', lw=.7)
            style(ax)
            numerator = lead_cell['n_called'] if j == 0 else lead_cell['copies'][1]
            denominator = data['totals'][group] if j == 0 else lead_cell['chromosomes']
            result = f'{numerator}/{denominator} ({numerator/denominator:.1%})' if denominator else 'NA (no usable calls)'
            ax.text(1, 1.04, 'Lead: '+result, transform=ax.transAxes, ha='right', fontsize=8.5, color='#40566d')
            if j == 1:
                missing = sum(cell['chromosomes'] == 0 for cell in cells)
                ax.text(.01, .06, f'{missing} sites with no calls: no frequency estimate', transform=ax.transAxes,
                        fontsize=8, color='#526579', bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': .9})
        axes[i, 0].set_ylabel(f'{display_group(group)}\nn = {data["totals"][group]}\nPercent', color=color, fontsize=10)
    for ax in axes[-1]:
        ax.set_xlabel(f'Distance from GWAS lead (kb; {summary["genome_build_user_declared"]})', labelpad=9)
    note = ('Left: usable participants / all participants in the group. Right: ALT copies / (2 × usable participants).\n'
            'Outlined diamond and vertical line identify the lead. ALT is the VCF alternate allele, not necessarily the GWAS effect or risk allele.\n')
    if summary.get('scan_counts', {}).get('high_quality_homozygous_reference_genotypes') == 0:
        note += 'No high-quality 0/0 calls were recorded in these windows; the observed fractions cannot be read as unbiased population frequencies.\n'
    note += 'These panels show recorded calls, not shared haplotypes, founder chromosomes, or tested group differences.'
    fig.text(.055, .118, note, fontsize=8.8, color='#526579', va='top', linespacing=1.65)
    caption = (prefix + f"{lead['locus_id']} was selected by highest pooled completeness of a unique exact biallelic SNP lead, "
               'with ties broken by GWAS lead-table order, independently of group frequency contrasts. '
               f'All {len(variants)} biallelic SNP records in this window are displayed; '
               f'{data["regional_excluded_records"]} other records are omitted. '
               'Each row is a participant group, not a chromosome copy. Positions use the actual genomic distance from the lead. '
               'Left: n_called / n_group. Right: ALT allele copies / called chromosomes, conditional on observed usable genotypes. '
               'Zero calls yield zero completeness but no allele-frequency estimate. Points are not connected and do not imply phase. '
               'The lead is outlined; ALT is not automatically the effect or risk allele. '
               'Missing and reference-call ascertainment can bias the fractions; no group association is inferred. '
               f'Layout inspiration: {PAPER_FIGURE}.')
    return save(plt, fig, out, 'regional_marker_comparison', prefix + title, caption,
                'synthetic_demo' if synthetic else 'cohort_run')
