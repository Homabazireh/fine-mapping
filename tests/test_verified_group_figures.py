"""Keep verified ID aliases and the plotted participant denominators consistent."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'code'))
import make_figures as figures


class VerifiedGroupFiguresTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.metadata = self.root / 'metadata.tsv'
        self.metadata.write_text('name\tDiabetes_Group\nsourceLR\tT2D\nexact\tControl\n')
        self.aliases = self.root / 'aliases.tsv'
        self.aliases.write_text('metadata_sample\tvcf_sample\tstatus\nsourceLR\ttarget\tconfirmed\n')
        table = self.root / 'qc.tsv'
        table.write_text('locus_id\tgroup\tn_group\thigh_quality_called_samples\n'
                         'L1\tT2D\t1\t1\nL1\tControl\t1\t0\n')
        self.summary = {
            'parameters': {'min_dp': 10, 'min_gq': 20},
            'inputs': {key: {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
                       for key, path in [('metadata', self.metadata), ('id_map', self.aliases)]},
            'group_qc': {
                'scope': 'exact_and_verified_id_matches', 'table': str(table),
                'source_category_counts': {'T2D': 1, 'Control': 1},
                'matched_category_counts': {'T2D': 1, 'Control': 1},
                'exact_match_category_counts': {'T2D': 0, 'Control': 1},
                'verified_alias_category_counts': {'T2D': 1, 'Control': 0},
                'matched_sample_count': 2, 'exact_match_sample_count': 1,
                'verified_alias_sample_count': 1, 'unmatched_vcf_count': 0,
            },
        }
        self.path = self.root / 'summary.json'

    def load(self):
        self.path.write_text(json.dumps(self.summary))
        return figures._load_preflight_qc(self.path, self.metadata)

    def test_verified_alias_denominators_and_provenance(self):
        _, _, rows = self.load()
        self.assertEqual([row['n_group'] for row in rows], ['1', '1'])
        self.aliases.write_text(self.aliases.read_text().replace('confirmed', 'pending'))
        with self.assertRaisesRegex(ValueError, 'ID-map SHA256'):
            self.load()

    def test_inconsistent_alias_counts_are_rejected(self):
        self.summary['group_qc']['verified_alias_category_counts']['T2D'] = 0
        with self.assertRaisesRegex(ValueError, 'exact and alias counts'):
            self.load()

    @unittest.skipUnless(importlib.util.find_spec('matplotlib'), 'Matplotlib is optional')
    def test_full_cohort_caption_counts_both_match_types(self):
        summary, _, rows = self.load()
        data = figures.prepare_genotype_qc_matrix(rows, ['L1'], ['T2D', 'Control'])
        output = figures._genotype_qc_figure(figures._pyplot(), data, summary, 1, self.root)
        self.assertIn('1 exact IDs and 1 explicitly verified aliases', output['caption'])
        self.assertIn('Diabetes: 1, Control: 1', output['caption'])
        self.assertNotIn('are not the plotting denominators', output['caption'])


if __name__ == '__main__':
    unittest.main()
