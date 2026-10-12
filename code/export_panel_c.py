#!/usr/bin/env python3
"""Export a standalone, high-resolution lead-completeness panel and report."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path

from build_report import select_run
from comparison_figures import load_comparison_data, PAPER_FIGURE
from group_labels import display_group
from make_figures import _load_run, _pyplot, _style_axis
from project_config import DEFAULT_CONFIG, load_config


def draw_panel_c(plt, data, summary):
    """Display group and pooled completeness at every selected GWAS lead."""
    from matplotlib.lines import Line2D

    colors = {'T2D': '#d55e00', 'Prediabetes': '#009e73', 'Control': '#0072b2', 'NORM': '#0072b2'}
    markers = ['o', '^', 's', 'v', 'P', 'X']
    groups, rows = data['groups'], data['coverage']
    figure, ax = plt.subplots(figsize=(13, max(8, 3.6 + len(rows)*.34)))
    figure.subplots_adjust(left=.255, right=.97, top=.77, bottom=.17)
    cutoff = summary['parameters']['min_call_rate']
    finite = [r for r in rows if math.isfinite(r['pooled'])]
    best = max(finite, key=lambda r: r['pooled']) if finite else None
    passed = sum(r['pooled'] >= cutoff for r in finite)
    figure.text(.045, .954, 'c   Usable lead genotypes across selected regions',
                fontsize=20, weight='bold', color='#172f46')
    figure.text(.045, .912, f"{data['total']} participants · {len(rows)} GWAS-selected regions · {summary['genome_build_user_declared']}",
                fontsize=12, color='#526579')
    params = summary['parameters']
    filters = ' · '.join(f'{label} ≥ {params[key]:g}' if params[key] else f'{label} filter disabled'
                        for key, label in [('min_dp', 'DP'), ('min_gq', 'GQ')])
    figure.text(.045, .877, f'Usable genotype: {filters}. Dashed line: {cutoff:.0%} pooled completeness requirement.',
                fontsize=10.5, color='#526579')
    offsets = [(i-(len(groups)-1)/2)*.19 for i in range(len(groups))]
    for j, row in enumerate(rows):
        ax.axhspan(j-.48, j+.48, color='#f5f7fa' if j % 2 == 0 else 'white', zorder=0)
        if row is best:
            ax.axhspan(j-.48, j+.48, color='#fff0ca', zorder=1)
    for i, group in enumerate(groups):
        ax.scatter([100*r['fractions'][i] for r in rows], [j+offsets[i] for j in range(len(rows))],
                    marker=markers[i % len(markers)], s=40, color=colors.get(group, '#64748b'),
                    edgecolor='white', linewidth=.5, zorder=4)
    ax.scatter([100*r['pooled'] for r in rows], range(len(rows)), marker='D', s=31,
                color='#233445', edgecolor='white', linewidth=.45, zorder=5)
    ax.axvline(100*cutoff, color='#a54536', lw=1.3, ls=(0, (4, 3)), zorder=2)
    labels = [f"{r['lead']['locus_id']}   {r['lead']['chrom']}:{int(r['lead']['lead_pos']):,}" for r in rows]
    ax.set_yticks(range(len(rows)), labels)
    ax.set_xticks([0, 20, 40, 60, 80, 100], ['0%', '20%', '40%', '60%', '80%', '100%'])
    ax.set_xlim(-2, 118)
    ax.set_ylim(len(rows)-.5, -.5)
    ax.set_xlabel('Participants with usable genotypes at the lead variant', fontsize=11, labelpad=11)
    ax.grid(axis='x', color='#dce3ea', linewidth=.65)
    _style_axis(ax)
    ax.tick_params(axis='y', length=0, labelsize=9.5, pad=9)
    ax.tick_params(axis='x', labelsize=10)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_bounds(0, 100)
    ax.text(109, 1.018, 'Pooled calls', transform=ax.get_xaxis_transform(), ha='center', fontsize=9, color='#526579')
    for j, row in enumerate(rows):
        value = f"{row['n_called']}/{data['total']}" if row['n_called'] is not None else 'NA'
        ax.text(109, j, value, ha='center', va='center', fontsize=9, weight='bold' if row is best else 'normal', color='#233445')
        if row is best:
            ax.get_yticklabels()[j].set_weight('bold')
    legend = [Line2D([], [], color=colors.get(g, '#64748b'), marker=markers[i % len(markers)], ls='',
                     markersize=7, label=f'{display_group(g)} (n = {data["totals"][g]})') for i, g in enumerate(groups)]
    legend.append(Line2D([], [], color='#233445', marker='D', ls='', markersize=6, label='All participants'))
    figure.legend(handles=legend, loc='upper center', bbox_to_anchor=(.58, .845), ncol=4,
                   frameon=False, fontsize=10, columnspacing=1.6, handletextpad=.5)
    if best:
        note = (f"Highlighted: {best['lead']['locus_id']} has the highest pooled coverage, "
                f"{best['n_called']}/{data['total']} ({best['pooled']:.1%}). "
                f'{passed}/{len(rows)} leads meet the pooled requirement.')
    else:
        note = 'No unique exact lead genotypes were available for pooled completeness estimates.'
    figure.text(.045, .09, note, fontsize=10, color='#40566d')
    figure.text(.045, .059, 'Dots measure genotype availability. They do not show haplotype frequencies or a tested diabetes association.',
                fontsize=9.5, color='#526579')
    figure.text(.045, .025, 'Real cohort · saved main-analysis tables · layout inspired by Figueroa et al., Nature Genetics (2024), Fig. 1',
                fontsize=8.5, color='#63758a')
    return figure


def export_panel_c(run_dir, out_dir, config=DEFAULT_CONFIG):
    selected, synthetic = select_run(load_config(config), run_dir=run_dir)
    if synthetic:
        raise ValueError('This export is for the real cohort; a synthetic run was selected')
    summary, leads, _, _ = _load_run(selected)
    data = load_comparison_data(selected, summary, leads)
    if data is None or not data['coverage']:
        raise ValueError('A completed run with regional genotype counts and selected leads is required')
    if not any(math.isfinite(row['pooled']) for row in data['coverage']):
        raise ValueError('No unique exact lead is available for a completeness estimate')
    out = Path(out_dir).resolve()
    if out.exists() and any(out.iterdir()):
        raise ValueError(f'Use a new or empty output folder: {out}')
    out.mkdir(parents=True, exist_ok=True)
    plt = _pyplot()
    with plt.rc_context({'font.family': 'DejaVu Sans', 'pdf.fonttype': 42, 'svg.fonttype': 'none'}):
        fig = draw_panel_c(plt, data, summary)
        for suffix in ['png', 'pdf', 'svg']:
            fig.savefig(out / f'panel_c.{suffix}', dpi=600, facecolor='white')
        plt.close(fig)
    svg_path = out / 'panel_c.svg'
    svg_path.write_text('\n'.join(line.rstrip() for line in svg_path.read_text().splitlines())+'\n')
    fields = ['locus_id', 'chrom', 'lead_pos', 'status', 'group', 'n_called', 'n_group', 'completeness']
    with (out / 'panel_c_data.tsv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for row in data['coverage']:
            for group in data['groups'] + ['ALL']:
                called = (row['n_called'] if group == 'ALL' else row['variant']['groups'][group]['n_called']) if row['variant'] else None
                total = data['total'] if group == 'ALL' else data['totals'][group]
                writer.writerow({**{k: row['lead'][k] for k in ['locus_id', 'chrom', 'lead_pos', 'status']},
                                 'group': group, 'n_called': called if called is not None else 'NA', 'n_group': total,
                                 'completeness': called/total if called is not None and total else 'NA'})
    finite = [r for r in data['coverage'] if math.isfinite(r['pooled'])]
    best = max(finite, key=lambda r: r['pooled']) if finite else None
    cutoff = summary['parameters']['min_call_rate']
    manifest = {
        'schema_version': 1, 'created_at_utc': datetime.now(timezone.utc).isoformat(), 'synthetic': False,
        'title': 'Panel c: usable lead genotypes across selected regions',
        'source_run': str(selected), 'source_run_completed_at_utc': summary['completed_at_utc'],
        'genome_build': summary['genome_build_user_declared'],
        'files': {'png': 'panel_c.png', 'pdf': 'panel_c.pdf', 'svg': 'panel_c.svg', 'data': 'panel_c_data.tsv'},
        'png_dpi': 600, 'group_counts': data['totals'], 'n_samples': data['total'], 'n_leads': len(leads),
        'min_call_rate': cutoff, 'lead_sites_passing': sum(r['pooled'] >= cutoff for r in finite),
        'max_high_quality_called': best['n_called'] if best else None,
        'max_high_quality_call_rate': best['pooled'] if best else None,
        'best_locus_id': best['lead']['locus_id'] if best else None,
        'best_lead': {k: best['variant'][k] for k in ['chrom', 'pos', 'ref', 'alt']} if best else None,
        'parameters': {k: summary['parameters'][k] for k in ['min_dp', 'min_gq', 'min_call_rate']},
        'locus_status_counts': summary['locus_status_counts'],
        'haplotype_carrier_tests_performed': summary['haplotype_carrier_tests_performed'],
        'high_quality_homozygous_reference_genotypes': summary['scan_counts']['high_quality_homozygous_reference_genotypes'],
        'preflight_disabled': True, 'source_paper_url': PAPER_FIGURE,
        'interpretation': 'Genotype completeness only; not haplotype frequencies, ancestry or group association results.',
    }
    manifest['sha256'] = {key: hashlib.sha256((out/name).read_bytes()).hexdigest() for key, name in manifest['files'].items()}
    (out / 'panel_c_manifest.json').write_text(json.dumps(manifest, indent=2, allow_nan=False)+'\n')
    from panel_c_report import write_report
    write_report(out, manifest)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', type=Path, required=True, help='Completed real cohort run')
    parser.add_argument('--out', type=Path, required=True, help='New or empty output folder')
    parser.add_argument('--config', default=str(DEFAULT_CONFIG))
    args = parser.parse_args()
    export_panel_c(args.run_dir, args.out, args.config)
    print(f'Panel c: 600-dpi PNG, vector PDF/SVG and focused report saved to {args.out.resolve()}')


if __name__ == '__main__':
    main()
