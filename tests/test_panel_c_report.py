"""Portable reporting of aggregate panel C data without loading cohort inputs."""
from __future__ import annotations

import base64
from html.parser import HTMLParser
import json
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'code'))
from panel_c_report import write_report


class _Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.images, self.downloads = [], []

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == 'img':
            self.images.append(attrs['src'])
        elif tag == 'a' and 'download' in attrs:
            self.downloads.append(attrs['href'])


class PanelCReportTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='panel_c_report_')
        self.addCleanup(temporary.cleanup)
        self.out = Path(temporary.name) / 'portable report'
        self.out.mkdir()
        self.png = b'\x89PNG\r\n\x1a\nopaque test image bytes'
        (self.out / 'panel_c.png').write_bytes(self.png)
        self.manifest = {
            'files': {'png': 'panel_c.png', 'pdf': 'panel_c.pdf',
                      'svg': 'panel_c.svg', 'data': 'panel_c_data.tsv'},
            'group_counts': {'Control': 50, 'Prediabetes': 48, 'T2D': 49},
            'n_leads': 20, 'n_samples': 147, 'min_call_rate': .9,
            'lead_sites_passing': 0, 'max_high_quality_called': 126,
            'max_high_quality_call_rate': 126 / 147,
            'parameters': {'min_dp': 10., 'min_gq': 20., 'min_call_rate': .9},
            'locus_status_counts': {'lead_below_min_call_rate': 19,
                                    'palindromic_lead_orientation_unconfirmed': 1},
            'haplotype_carrier_tests_performed': 0,
            'high_quality_homozygous_reference_genotypes': 0,
            'preflight_disabled': True, 'synthetic': False,
            'source_paper_url': 'https://www.nature.com/articles/s41588-024-01719-5/figures/1',
            'source_run_completed_at_utc': '2026-10-12T17:33:17+00:00', 'genome_build': 'GRCh38',
            'private_source_path': '/private/participant_data.tsv',
        }

    def test_one_embedded_image_and_relative_downloads_remain_portable(self):
        for name in ['panel_c.pdf', 'panel_c.svg', 'panel_c_data.tsv']:
            (self.out / name).write_text('Existing exported asset')
        (self.out / 'panel_c_manifest.json').write_text(json.dumps(self.manifest))
        original_files = {path.name: path.read_bytes() for path in self.out.iterdir()}
        result = write_report(self.out, self.manifest)
        self.assertEqual(result, self.out / 'report.html')
        html = result.read_text()
        parser = _Links()
        parser.feed(html)
        self.assertEqual(parser.images, ['data:image/png;base64,' + base64.b64encode(self.png).decode()])
        self.assertEqual(set(parser.downloads), {*self.manifest['files'].values(),
                                               'panel_c_manifest.json', 'REPORT.txt'})
        self.assertIn('None of the 20 lead sites meets 90% completeness', html)
        self.assertIn('49 Diabetes, 48 Prediabetes, 50 Control', html)
        self.assertIn('126/147', html)
        self.assertIn('A genotype is a DNA call at one position.', html)
        self.assertIn('read-depth (DP) and call-quality (GQ)', html)
        self.assertIn('Labels L0001–L0020 identify regions ranked by GWAS significance', html)
        self.assertIn('0 usable haplotype panels', html)
        self.assertIn('0 haplotype association tests', html)
        self.assertIn('The separate preliminary VCF check was skipped.', html)
        self.assertIn('retained its internal genotype and phase requirements', html)
        self.assertNotIn(self.manifest['private_source_path'], html)
        self.assertNotIn(str(self.out), html)
        relocated = self.out.parent / 'downloaded elsewhere'
        self.out.rename(relocated)
        for href in parser.downloads:
            self.assertTrue((relocated / href).is_file())
        for name, content in original_files.items():
            self.assertEqual((relocated / name).read_bytes(), content)
        self.assertIn('Missing genotypes remain missing.', (relocated / 'REPORT.txt').read_text())

    def test_inconsistent_or_synthetic_counts_do_not_produce_a_real_cohort_report(self):
        for change in [{'synthetic': True}, {'n_samples': 146},
                       {'max_high_quality_call_rate': .5}, {'min_call_rate': .8}]:
            with self.subTest(change=change):
                with self.assertRaises(ValueError):
                    write_report(self.out, {**self.manifest, **change})
                self.assertFalse((self.out / 'report.html').exists())


if __name__ == '__main__':
    unittest.main()
