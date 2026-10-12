"""Write a portable, single-figure explanation of recorded lead-genotype completeness."""
from __future__ import annotations

import base64
from html import escape
from pathlib import Path
from urllib.parse import urlsplit

from group_labels import display_group, ordered_groups


def write_report(out: Path, manifest: dict) -> Path:
    """Embed the exported PNG and write aggregate-only HTML and plain-text reports.

    The exporter supplies counts from its completed-run tables. This function
    neither reads a VCF nor loads a preflight summary. Relative download links
    keep the report usable when its containing artifact directory is copied.
    """
    out = Path(out)
    if manifest.get('synthetic') is not False:
        raise ValueError('The panel C cohort report requires a non-synthetic manifest.')
    files = manifest['files']
    for key in ('png', 'pdf', 'svg', 'data'):
        filename = files[key]
        if not isinstance(filename, str) or Path(filename).name != filename:
            raise ValueError('Panel C report downloads must be filenames in the output directory.')
    png = (out / files['png']).read_bytes()
    if not png.startswith(b'\x89PNG\r\n\x1a\n'):
        raise ValueError('The panel C image is not a PNG file.')

    groups = manifest['group_counts']
    n_samples = manifest['n_samples']
    if n_samples <= 0 or sum(groups.values()) != n_samples:
        raise ValueError('Participant group counts must sum to the positive sample count.')
    n_leads, passing = manifest['n_leads'], manifest['lead_sites_passing']
    best = manifest['max_high_quality_called']
    if not 0 <= passing <= n_leads or not 0 <= best <= n_samples:
        raise ValueError('Lead-site summary counts are inconsistent.')
    parameters = manifest['parameters']
    cutoff = manifest['min_call_rate']
    if cutoff != parameters['min_call_rate'] or not 0 <= cutoff <= 1:
        raise ValueError('Genotype-completeness thresholds are inconsistent.')
    if abs(manifest['max_high_quality_call_rate'] - best / n_samples) > 1e-7:
        raise ValueError('The best lead call rate disagrees with its count.')

    paper_url = manifest['source_paper_url']
    if urlsplit(paper_url).scheme not in {'https', 'http'}:
        raise ValueError('The paper citation must be an HTTP(S) URL.')
    cutoff_text = f'{cutoff:.0%}'
    title = (f'None of the {n_leads} lead sites meets {cutoff_text} completeness'
             if passing == 0 else f'{passing} of {n_leads} lead sites meet {cutoff_text} completeness')
    group_text = ', '.join(f'{groups[group]} {display_group(group)}' for group in ordered_groups(groups))
    quality = (f'DP ≥ {parameters["min_dp"]:g} and GQ ≥ {parameters["min_gq"]:g}')
    result = (f'This completed analysis included all {n_samples} participants: {group_text}. '
              f'At the best-covered lead, {best}/{n_samples} participants '
              f'({best / n_samples:.1%}) had a genotype meeting {quality}.')
    guide = ('A genotype is a DNA call at one position. A usable genotype meets the recorded '
             'read-depth (DP) and call-quality (GQ) requirements. '
             f'Labels L0001–L{n_leads:04d} identify regions ranked by GWAS significance; '
             'they are not participant IDs. Each row represents one GWAS lead site. '
             'Colored dots show the share of participants '
             'with usable genotypes within each group. Black diamonds combine all participants; '
             f'the counts on the right use {n_samples} as their denominator. The dashed line '
             f'marks the {cutoff_text} completeness requirement. Positions further right mean '
             'more complete genotype data, not greater diabetes risk or a stronger genetic effect.')
    states = manifest['locus_status_counts']
    analysed = states.get('analysed', 0)
    tests = manifest['haplotype_carrier_tests_performed']
    outcome = (f'The completed run produced {analysed} usable haplotype panels and performed '
               f'{tests} haplotype association tests.')
    if analysed == 0:
        outcome += (' No haplotype frequencies or differences in those frequencies between '
                    'groups could be estimated. '
                    'This does not establish that the groups have identical haplotypes.')
    skipped = []
    for key, label in (
        ('lead_below_min_call_rate', 'insufficient usable lead genotypes'),
        ('palindromic_lead_orientation_unconfirmed', 'unconfirmed palindromic allele orientation'),
    ):
        if states.get(key):
            skipped.append(f'{states[key]} for {label}')
    if skipped:
        outcome += ' Recorded region exclusions: ' + '; '.join(skipped) + '.'
    source = ('The separate preliminary VCF check was skipped. '
              if manifest['preflight_disabled'] else '')
    source += ('This figure uses tables already saved by the completed analysis; creating it '
               'does not rescan the VCF. The analysis retained its internal genotype and phase '
               'requirements. Missing genotypes remain missing.')
    reference_note = ''
    if manifest['high_quality_homozygous_reference_genotypes'] == 0:
        reference_note = ('No high-quality homozygous-reference (0/0) calls were recorded in '
                          'the selected windows. The figure describes the available calls; '
                          'it cannot identify the cause of this absence.')
    citation = ('Visual organization inspired by Figueroa et al. (2024), Figure 1. '
                'Every plotted value comes from this cohort; the figure shows genotype '
                'completeness, not the paper\'s shared chromosome patterns.')
    downloads = [
        (files['png'], 'PNG · 600 dpi'), (files['pdf'], 'PDF · scalable'),
        (files['svg'], 'SVG · editable'), (files['data'], 'Plotted data · TSV'),
        ('panel_c_manifest.json', 'Figure provenance · JSON'), ('REPORT.txt', 'Plain-text report'),
    ]
    links = '\n'.join(f'<a href="{escape(path, quote=True)}" download>{escape(label)}</a>'
                      for path, label in downloads)
    rows = ''.join(f'<tr><td>{escape(display_group(group))}</td><td>{groups[group]}</td></tr>'
                   for group in ordered_groups(groups))
    note_html = f'<p>{escape(reference_note)}</p>' if reference_note else ''
    image_data = base64.b64encode(png).decode('ascii')
    content = f'''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{escape(title)}</title>
<style>
:root {{ color-scheme: light; font-family: system-ui, -apple-system, sans-serif; color: #203347; background: #f5f7fa; }}
body {{ max-width: 1100px; margin: 0 auto; padding: 40px 24px 56px; line-height: 1.6; }}
main {{ background: white; padding: 32px; border-radius: 12px; box-shadow: 0 2px 12px #2033470d; }}
h1 {{ font-size: clamp(1.6rem, 3vw, 2.35rem); line-height: 1.2; margin: 0 0 20px; }}
h2 {{ font-size: 1.25rem; margin-top: 30px; }}
.kicker {{ color: #5d6e80; font-size: .85rem; letter-spacing: .08em; text-transform: uppercase; }}
figure {{ margin: 28px 0; }}
figure img {{ display: block; width: 100%; height: auto; }}
figcaption {{ font-size: .95rem; color: #425569; margin-top: 16px; }}
table {{ border-collapse: collapse; min-width: 280px; }}
th, td {{ padding: 8px 18px 8px 0; border-bottom: 1px solid #dce3eb; text-align: left; }}
.downloads {{ display: flex; flex-wrap: wrap; gap: 10px; margin: 20px 0; }}
a {{ color: #245b89; }}
.downloads a {{ border: 1px solid #cbd6e2; border-radius: 6px; padding: 7px 12px; text-decoration: none; }}
.provenance {{ color: #526579; font-size: .9rem; border-top: 1px solid #dce3eb; margin-top: 30px; padding-top: 18px; }}
@media (max-width: 600px) {{ body {{ padding: 12px; }} main {{ padding: 18px; }} }}
@media print {{ body, main {{ background: white; padding: 0; box-shadow: none; }} .downloads {{ display: none; }} }}
</style>
</head>
<body><main>
<p class="kicker">Panel C · Lead-genotype completeness</p>
<h1>{escape(title)}</h1>
<p>{escape(result)}</p>
<table aria-label="Participants included in the analysis"><thead><tr><th>Group</th><th>Participants</th></tr></thead><tbody>{rows}</tbody></table>
<figure>
<img src="data:image/png;base64,{image_data}" alt="Lead-genotype completeness for Diabetes, Prediabetes and Control, with pooled counts and the completeness threshold">
<figcaption>{escape(guide)}</figcaption>
</figure>
<nav class="downloads" aria-label="Download figure and report files">{links}</nav>
<h2>What the result means</h2>
<p>{escape(outcome)}</p>
<p>{escape(source)}</p>
{note_html}
<div class="provenance">
<p>{escape(citation)} <a href="{escape(paper_url, quote=True)}">View the source figure</a>.</p>
<p>Source run completed (UTC): {escape(str(manifest['source_run_completed_at_utc']))}. Genome build: {escape(str(manifest['genome_build']))}.</p>
</div>
</main></body></html>
'''
    out.mkdir(parents=True, exist_ok=True)
    report = out / 'report.html'
    report.write_text(content, encoding='utf-8')
    paragraphs = [title.upper(), result, 'HOW TO READ PANEL C', guide,
                  'WHAT THE RESULT MEANS', outcome, source, reference_note,
                  'SOURCE', citation, paper_url,
                  f'Source run completed (UTC): {manifest["source_run_completed_at_utc"]}. '
                  f'Genome build: {manifest["genome_build"]}.',
                  'FILES\n' + '\n'.join(path + ': ' + label for path, label in downloads)]
    (out / 'REPORT.txt').write_text('\n\n'.join(p for p in paragraphs if p) + '\n', encoding='utf-8')
    return report
