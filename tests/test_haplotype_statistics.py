"""Participant-level inference, provenance gates, and independent exact-test checks."""
from __future__ import annotations

from collections import Counter
import contextlib
import csv
import io
from itertools import combinations
import math
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "code"))

import compare_haplotypes as pipeline
from make_demo import make_demo

try:
    from scipy.stats import fisher_exact
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False


class ExactCarrierTests(unittest.TestCase):
    def test_overall_matches_exhaustive_participant_label_assignments(self):
        # Enumerate all assignments of four carriers to nine distinct people.
        # This independently checks the table-weight computation and its ties.
        observed = (0, 1, 3)
        outcomes = Counter()
        for carriers in combinations(range(9), 4):
            outcomes[tuple(sum(person // 3 == group for person in carriers) for group in range(3))] += 1
        expected = sum(count for count in outcomes.values() if count <= outcomes[observed]) / math.comb(9, 4)
        actual, n_tables = pipeline.overall_carrier_test(list(observed), [3, 3, 3])
        self.assertEqual(actual, expected)
        self.assertEqual(n_tables, len(outcomes))

    def test_three_group_exact_matches_R_fisher_reference(self):
        # Reference: stats::fisher.test(rbind(c(3,7),c(8,2),c(6,4))) in R.
        actual, n_tables = pipeline.overall_carrier_test([3, 8, 6], [10, 10, 10])
        self.assertAlmostEqual(actual, 0.10299286447002098, places=13)
        self.assertGreater(n_tables, 1)
        self.assertEqual(pipeline.overall_carrier_test([0, 1, 2], [2, 2, 2])[0], 0.6)

    def test_overall_is_invariant_to_group_order_and_carrier_complement(self):
        reference = pipeline.overall_carrier_test([2, 7, 3], [8, 10, 5])
        self.assertEqual(reference, pipeline.overall_carrier_test([3, 2, 7], [5, 8, 10]))
        self.assertEqual(reference, pipeline.overall_carrier_test([6, 3, 2], [8, 10, 5]))
        self.assertEqual(pipeline.overall_carrier_test([0, 0, 0], [8, 10, 5]), (1.0, 1))

    @unittest.skipUnless(SCIPY_AVAILABLE, "SciPy is optional")
    def test_two_group_enumeration_matches_independent_scipy_fisher(self):
        for a, n_a, b, n_b in [(1, 10, 11, 14), (0, 10, 7, 12), (3, 4, 1, 4), (49, 49, 48, 48)]:
            with self.subTest(counts=(a, n_a, b, n_b)):
                expected = fisher_exact([[a, n_a - a], [b, n_b - b]]).pvalue
                self.assertAlmostEqual(pipeline.overall_carrier_test([a, b], [n_a, n_b])[0], expected, places=13)

    @unittest.skipUnless(SCIPY_AVAILABLE, "SciPy is optional")
    def test_pairwise_direction_and_exact_intervals_match_R_fisher(self):
        odds, low, high, p = pipeline.carrier_test(1, 10, 11, 14)
        # R fisher.test uses less stringent numerical root tolerances than SciPy.
        self.assertAlmostEqual(p, 0.0027594561852200836, places=13)
        self.assertAlmostEqual(odds, 0.037233117250979164, delta=5e-5)
        self.assertAlmostEqual(low, 0.000643828388365295, delta=1e-5)
        self.assertAlmostEqual(high, 0.4258840381150771, delta=1e-4)
        reverse = pipeline.carrier_test(11, 14, 1, 10)
        self.assertAlmostEqual(odds * reverse[0], 1.0, places=9)
        self.assertAlmostEqual(low * reverse[2], 1.0, places=9)
        self.assertAlmostEqual(high * reverse[1], 1.0, places=9)
        self.assertAlmostEqual(p, reverse[3], places=13)

    @unittest.skipUnless(SCIPY_AVAILABLE, "SciPy is optional")
    def test_zero_cells_preserve_boundary_odds_and_finite_exact_bound(self):
        odds, low, high, p = pipeline.carrier_test(0, 10, 7, 12)
        self.assertEqual((odds, low), (0.0, 0.0))
        self.assertAlmostEqual(high, 0.5500622061965632, delta=1e-4)
        self.assertAlmostEqual(p, 0.005347593582887695, places=13)

    def test_score_difference_interval_handles_boundaries_and_reversal(self):
        for counts in [(0, 10, 0, 10), (0, 10, 10, 10), (8, 10, 3, 10), (49, 49, 48, 48)]:
            difference, low, high = pipeline.carrier_difference_interval(*counts)
            a, n_a, b, n_b = counts
            reverse = pipeline.carrier_difference_interval(b, n_b, a, n_a)
            self.assertAlmostEqual(difference, a / n_a - b / n_b)
            self.assertLess(low, high)
            self.assertLessEqual(low, difference)
            self.assertGreaterEqual(high, difference)
            self.assertGreaterEqual(low, -1)
            self.assertLessEqual(high, 1)
            self.assertAlmostEqual(difference, -reverse[0])
            self.assertAlmostEqual(low, -reverse[2])
            self.assertAlmostEqual(high, -reverse[1])

    def test_bh_family_does_not_include_untested_rows(self):
        self.assertEqual(pipeline.bh_adjust([0.01, 0.02, None, 0.5]), [0.03, 0.03, None, 0.5])
        self.assertEqual(pipeline.bh_adjust([None, None]), [None, None])


class StatisticalWorkflowTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="test_haplotype_statistics_")
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)
        self.fixtures = make_demo(self.folder / "inputs")
        self.out = self.folder / "results"
        self.args = pipeline.argument_parser().parse_args([
            "--vcf", str(self.fixtures / "cohort.vcf"),
            "--gwas", str(self.fixtures / "gwas.tsv"),
            "--samples", str(self.fixtures / "samples.tsv"),
            "--out", str(self.out), "--genome-build", "GRCh38",
            "--test-associations", "--min-test-samples", "1", "--min-carriers", "1",
            "--min-resolved-fraction", "0.5",
        ])

    def run_analysis(self):
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            return pipeline.run(self.args)

    def read_output(self, filename):
        with (self.out / filename).open() as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_requested_tests_without_reference_provenance_remain_unavailable(self):
        with mock.patch.object(pipeline, "carrier_test", side_effect=AssertionError("Invalid inference")), \
             mock.patch.object(pipeline, "overall_carrier_test", side_effect=AssertionError("Invalid inference")):
            summary = self.run_analysis()
        overall = self.read_output("haplotype_overall_comparisons.tsv")
        pairwise = self.read_output("haplotype_comparisons.tsv")
        self.assertTrue(overall and pairwise)
        self.assertTrue(summary["association_tests_requested"])
        self.assertFalse(summary["association_tests_permitted"])
        self.assertEqual(summary["haplotype_overall_tests_performed"], 0)
        self.assertEqual(summary["haplotype_pairwise_tests_performed"], 0)
        for row in overall + pairwise:
            self.assertEqual(row["test_status"], "not_tested_reference_provenance")
        for row in overall:
            self.assertEqual((row["p_fisher_overall"], row["q_BH_all_overall_tests"]), ("NA", "NA"))
        for row in pairwise:
            for field in ["p_fisher_two_sided", "conditional_odds_ratio_a_vs_b", "CI95_low", "CI95_high",
                          "carrier_frequency_difference_a_minus_b", "carrier_difference_CI95_low",
                          "carrier_difference_CI95_high", "q_BH_all_pairwise_tests"]:
                self.assertEqual(row[field], "NA")

    def test_failed_loci_produce_empty_tables_not_null_associations(self):
        self.args.min_dp = 1e9
        summary = self.run_analysis()
        self.assertEqual(summary["locus_status_counts"], {"lead_below_min_call_rate": 2})
        self.assertEqual(self.read_output("haplotype_overall_comparisons.tsv"), [])
        self.assertEqual(self.read_output("haplotype_comparisons.tsv"), [])
        self.assertEqual(summary["haplotype_overall_tests_performed"], 0)
        self.assertEqual(summary["haplotype_pairwise_tests_performed"], 0)

    @unittest.skipUnless(SCIPY_AVAILABLE, "SciPy is optional")
    def test_global_and_pairwise_families_use_separate_BH_corrections(self):
        self.args.reference_genotypes_confirmed = True  # Only generated, fully known calls.
        summary = self.run_analysis()
        overall = self.read_output("haplotype_overall_comparisons.tsv")
        pairwise = self.read_output("haplotype_comparisons.tsv")
        for rows, p_key, q_key in [(overall, "p_fisher_overall", "q_BH_all_overall_tests"),
                                   (pairwise, "p_fisher_two_sided", "q_BH_all_pairwise_tests")]:
            tested = sorted([row for row in rows if row[p_key] != "NA"], key=lambda row: float(row[p_key]))
            self.assertTrue(tested)
            # Independently evaluate the BH definition within this family.
            for rank, row in enumerate(tested, 1):
                expected = min([1.0] + [len(tested) * float(later[p_key]) / j
                                       for j, later in enumerate(tested[rank - 1:], rank)])
                self.assertAlmostEqual(float(row[q_key]), expected, places=7)
        for row in pairwise:
            self.assertEqual(row["q_BH_all_haplotype_tests"], row["q_BH_all_pairwise_tests"])
            if row["test_status"] == "tested_unadjusted":
                expected = float(row["carrier_frequency_a"]) - float(row["carrier_frequency_b"])
                self.assertAlmostEqual(float(row["carrier_frequency_difference_a_minus_b"]), expected, places=7)
        self.assertEqual(summary["haplotype_overall_tests_performed"],
                         sum(row["test_status"] == "tested_unadjusted" for row in overall))
        self.assertEqual(summary["haplotype_pairwise_tests_performed"],
                         sum(row["test_status"] == "tested_unadjusted" for row in pairwise))

    def test_resolution_and_carrier_count_gates_apply_before_inference(self):
        self.args.reference_genotypes_confirmed = True
        self.args.min_test_samples = 10
        self.assertEqual(pipeline.carrier_test_status(self.args, [1, 1, 1], [4, 4, 4], [4, 4, 4]),
                         "too_few_resolved_samples")
        self.args.min_test_samples = 1
        self.args.min_resolved_fraction = 0.8
        self.assertEqual(pipeline.carrier_test_status(self.args, [2, 2, 2], [4, 4, 4], [4, 4, 10]),
                         "low_resolved_fraction")
        self.assertEqual(pipeline.carrier_test_status(self.args, [4, 4, 4], [4, 4, 4], [4, 4, 4]),
                         "no_carrier_variation")
        self.args.min_carriers = 3
        self.assertEqual(pipeline.carrier_test_status(self.args, [1, 1, 0], [4, 4, 4], [4, 4, 4]),
                         "too_few_carriers_or_noncarriers")


if __name__ == "__main__":
    unittest.main()
