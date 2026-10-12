"""Scientific count and selection invariants for completed-run comparison plots."""
from __future__ import annotations

from copy import deepcopy
import math
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'code'))

from comparison_figures import prepare_comparison_data


class ComparisonDataTests(unittest.TestCase):
    def setUp(self):
        # Unequal group sizes expose accidental averaging of group percentages.
        self.summary = {'group_counts': {'T2D': 4, 'Prediabetes': 2, 'Control': 2}}

    @staticmethod
    def lead(ident, position):
        return {'locus_id': ident, 'chrom': 'chr1', 'lead_pos': str(position),
                'start': str(position - 50), 'end': str(position + 50)}

    def variant(self, locus, position, cells, *, alt='G', kind='SNP', distance=0, exact=True):
        """Describe explicit group counts as (called people, allele-copy counts)."""
        rows = []
        for group, (called, copies) in cells.items():
            for index, count in enumerate(copies):
                rows.append({
                    'locus_id': locus, 'variant_id': f'chr1:{position}:A:{alt}',
                    'chrom': 'chr1', 'pos': str(position), 'ref': 'A', 'alt': alt,
                    'type': kind, 'distance_to_lead_bp': str(distance),
                    'is_exact_lead_match': '1' if exact else '0', 'group': group,
                    'n_group': str(self.summary['group_counts'][group]),
                    'n_called': str(called), 'called_chromosomes': str(2 * called),
                    'allele_index': str(index), 'allele_copies': str(count),
                    # Deliberately misleading cached fractions must be ignored.
                    'allele_frequency': '0.987', 'call_rate': '0.987',
                })
        return rows

    def basic_rows(self):
        return self.variant('L1', 100, {
            'T2D': (2, [3, 1]), 'Prediabetes': (1, [1, 1]), 'Control': (0, [0, 0]),
        })

    def test_zero_calls_mean_zero_completeness_and_undefined_allele_frequency(self):
        data = prepare_comparison_data(self.basic_rows(), self.summary, [self.lead('L1', 100)])
        cells = data['example']['variant']['groups']
        self.assertEqual(cells['Control']['completeness'], 0)
        self.assertTrue(math.isnan(cells['Control']['alt_fraction']))
        self.assertEqual(cells['T2D']['completeness'], 0.5)
        self.assertEqual(cells['T2D']['alt_fraction'], 0.25)
        self.assertEqual(cells['Prediabetes']['alt_fraction'], 0.5)
        self.assertEqual(data['example']['n_called'], 3)
        self.assertEqual(data['example']['pooled'], 3 / 8)

    def test_example_uses_pooled_counts_then_lead_order_independently_of_phenotype_contrasts(self):
        leads = [self.lead('low', 100), self.lead('first_tie', 200), self.lead('second_tie', 300)]
        # The low-coverage lead has extreme group allele contrasts. The two
        # higher-coverage leads tie despite different group call distributions.
        rows = self.variant('low', 100, {
            'T2D': (1, [2, 0]), 'Prediabetes': (1, [0, 2]), 'Control': (1, [0, 2]),
        }) + self.variant('first_tie', 200, {
            'T2D': (4, [4, 4]), 'Prediabetes': (0, [0, 0]), 'Control': (0, [0, 0]),
        }) + self.variant('second_tie', 300, {
            'T2D': (0, [0, 0]), 'Prediabetes': (2, [0, 4]), 'Control': (2, [4, 0]),
        })
        for ordered_rows in (rows, list(reversed(rows))):
            data = prepare_comparison_data(ordered_rows, self.summary, leads)
            self.assertEqual(data['example']['lead']['locus_id'], 'first_tie')
            self.assertEqual(data['example']['pooled'], 0.5)
        reordered = prepare_comparison_data(rows, self.summary, [leads[2], leads[1], leads[0]])
        self.assertEqual(reordered['example']['lead']['locus_id'], 'second_tie')
        # Change every observed allele contrast without changing availability.
        no_contrasts = deepcopy(rows)
        for row in no_contrasts:
            row['allele_copies'] = row['n_called']
        unchanged = prepare_comparison_data(no_contrasts, self.summary, leads)
        self.assertEqual(unchanged['example']['lead']['locus_id'], 'first_tie')

    def test_called_chromosomes_must_be_twice_called_participants(self):
        rows = self.basic_rows()
        rows[0]['called_chromosomes'] = '6'
        with self.assertRaisesRegex(ValueError, 'denominators are inconsistent'):
            prepare_comparison_data(rows, self.summary, [self.lead('L1', 100)])

    def test_allele_counts_must_sum_to_the_call_denominator(self):
        rows = self.basic_rows()
        rows[0]['allele_copies'] = '2'  # 2 + 1 does not account for 4 chromosomes.
        with self.assertRaisesRegex(ValueError, 'copies do not sum'):
            prepare_comparison_data(rows, self.summary, [self.lead('L1', 100)])

    def test_duplicate_or_missing_allele_rows_are_rejected(self):
        for rows in [self.basic_rows() + [self.basic_rows()[0]], self.basic_rows()[1:]]:
            with self.subTest(rows=len(rows)):
                with self.assertRaises(ValueError):
                    prepare_comparison_data(rows, self.summary, [self.lead('L1', 100)])

    def test_alleles_of_one_group_must_agree_on_called_participants(self):
        rows = self.basic_rows()
        rows[1].update(n_called='1', called_chromosomes='2')
        with self.assertRaisesRegex(ValueError, 'inconsistent denominators'):
            prepare_comparison_data(rows, self.summary, [self.lead('L1', 100)])

    def test_every_variant_requires_all_groups_even_groups_without_calls(self):
        rows = [row for row in self.basic_rows() if row['group'] != 'Control']
        with self.assertRaisesRegex(ValueError, 'missing a participant group'):
            prepare_comparison_data(rows, self.summary, [self.lead('L1', 100)])

    def test_group_size_and_unknown_groups_or_loci_are_rejected(self):
        for change in [{'n_group': '5'}, {'group': 'unknown'}, {'locus_id': 'unknown'}]:
            with self.subTest(change=change):
                rows = self.basic_rows()
                rows[0].update(change)
                with self.assertRaises(ValueError):
                    prepare_comparison_data(rows, self.summary, [self.lead('L1', 100)])

    def test_multiallelic_lead_has_coverage_but_is_not_a_regional_example(self):
        leads = [self.lead('multi', 100), self.lead('biallelic', 200)]
        multi_cells = {'T2D': (4, [2, 3, 3]), 'Prediabetes': (2, [1, 1, 2]), 'Control': (2, [0, 2, 2])}
        rows = self.variant('multi', 100, multi_cells, alt='G,T')
        rows += self.variant('biallelic', 200, {
            'T2D': (2, [3, 1]), 'Prediabetes': (1, [1, 1]), 'Control': (0, [0, 0]),
        })
        rows += self.variant('biallelic', 210, multi_cells, alt='G,T', distance=10, exact=False)
        rows += self.variant('biallelic', 205, {
            'T2D': (1, [1, 1]), 'Prediabetes': (1, [1, 1]), 'Control': (1, [1, 1]),
        }, alt='AT', kind='indel', distance=5, exact=False)
        data = prepare_comparison_data(rows, self.summary, leads)
        self.assertEqual(data['coverage'][0]['n_called'], 8)
        self.assertEqual(data['coverage'][0]['pooled'], 1)
        self.assertTrue(math.isnan(data['coverage'][0]['variant']['groups']['T2D']['alt_fraction']))
        self.assertEqual(data['example']['lead']['locus_id'], 'biallelic')
        self.assertEqual([v['pos'] for v in data['regional']], [200])
        self.assertEqual(data['regional_total_records'], 3)
        self.assertEqual(data['regional_excluded_records'], 2)

    def test_ambiguous_exact_leads_have_no_coverage_estimate(self):
        cells = {'T2D': (2, [3, 1]), 'Prediabetes': (1, [1, 1]), 'Control': (0, [0, 0])}
        rows = self.variant('L1', 100, cells) + self.variant('L1', 100, cells, alt='T')
        data = prepare_comparison_data(rows, self.summary, [self.lead('L1', 100)])
        self.assertIsNone(data['coverage'][0]['n_called'])
        self.assertTrue(all(math.isnan(value) for value in data['coverage'][0]['fractions']))
        self.assertIsNone(data['example'])
        self.assertEqual(data['regional'], [])

    def test_all_uncalled_lead_is_reported_but_not_selected_for_regional_display(self):
        rows = self.variant('L1', 100, {g: (0, [0, 0]) for g in self.summary['group_counts']})
        data = prepare_comparison_data(rows, self.summary, [self.lead('L1', 100)])
        self.assertEqual(data['coverage'][0]['fractions'], [0, 0, 0])
        self.assertEqual(data['coverage'][0]['pooled'], 0)
        self.assertIsNone(data['example'])


if __name__ == '__main__':
    unittest.main()
