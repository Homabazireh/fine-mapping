"""Synthetic regressions for descriptive contrast pictures and unavailable cells."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "code"))
from plot_haplotype_contrasts import draw_contrasts, draw_region_contrasts


@unittest.skipUnless(importlib.util.find_spec("matplotlib"), "Matplotlib optional")
class ContrastRenderingTests(unittest.TestCase):
    def setUp(self):
        from matplotlib import pyplot as plt
        self.addCleanup(plt.close, "all")
        self.summary = {"synthetic": True, "selected_loci": 4,
                        "group_counts": {"T2D": 20, "Prediabetes": 20, "Control": 20},
                        "parameters": {"min_call_rate": .9}}

    def rows(self, n_loci=3, n_haps=5):
        rows = []
        for locus in range(1, n_loci + 1):
            for hap in range(1, n_haps + 1):
                for group in ["Diabetes", "Prediabetes"]:
                    rows.append(dict(locus_id=f"L{locus:04d}", chrom="chr1", lead_pos=100*locus,
                        haplotype_id=f"L{locus:04d}_H{hap:03d}", contrast=f"{group} vs Control",
                        comparison_status="descriptive_only", reason="not_tested_reference_provenance",
                        n_resolved_case=20, n_resolved_control=20, carriers_case=hap, carriers_control=1,
                        descriptive_carrier_frequency_difference=(hap-1)/20,
                        carrier_frequency_difference=None, CI95_low=None, CI95_high=None,
                        p_fisher_two_sided=None, q_BH_all_pairwise_tests=None))
        return rows

    def test_overview_distinguishes_observed_zero_missing_and_excluded_regions(self):
        rows = self.rows()
        rows[1]["comparison_status"] = "not_estimable"
        rows[1]["descriptive_carrier_frequency_difference"] = None
        rows += [dict(rows[0], locus_id="L0004", haplotype_id=None,
                      comparison_status="not_estimable", reason="lead_below_min_call_rate",
                      descriptive_carrier_frequency_difference=None),
                 dict(rows[1], locus_id="L0004", haplotype_id=None,
                      comparison_status="not_estimable", reason="lead_below_min_call_rate")]
        figure = draw_contrasts(self.summary, rows)
        self.assertEqual(len(figure.axes), 4)  # Three observed regions plus one shared colorbar.
        image = figure.axes[0].images[0].get_array()
        self.assertEqual(float(image[0, 0]), 0)
        self.assertFalse(bool(image.mask[0, 0]))
        self.assertTrue(bool(image.mask[0, 1]))
        labels = [t.get_text() for t in figure.axes[0].texts]
        self.assertIn("0.00", labels)
        self.assertIn("NA", labels)
        self.assertTrue(any("DESCRIPTIVE ONLY" in t.get_text() for t in figure.texts))
        self.assertTrue(any("1 additional regions" in t.get_text() for t in figure.texts))
        self.assertTrue(all(len(ax.containers) == 0 for ax in figure.axes))

    def test_regional_points_include_counts_without_invented_uncertainty(self):
        figure = draw_region_contrasts(self.summary, self.rows(n_loci=1, n_haps=2))
        self.assertEqual(len(figure.axes), 2)
        for ax in figure.axes:
            self.assertEqual(sum(line.get_marker() == "D" for line in ax.lines), 2)
            self.assertEqual(len(ax.containers), 0)
            self.assertTrue(any("2/20 vs 1/20" in t.get_text() and "+5.00 pp" in t.get_text() for t in ax.texts))
        self.assertTrue(any("DESCRIPTIVE ONLY" in t.get_text() for t in figure.texts))

    def test_tested_regional_point_keeps_its_interval_and_q(self):
        rows = self.rows(n_loci=1, n_haps=1)
        rows[0].update(comparison_status="tested", carrier_frequency_difference=.1,
                       CI95_low=-.05, CI95_high=.2, q_BH_all_pairwise_tests=.12)
        figure = draw_region_contrasts(self.summary, rows)
        self.assertEqual(len(figure.axes[0].containers), 1)
        self.assertEqual(len(figure.axes[1].containers), 0)
        self.assertTrue(any("q=0.12" in t.get_text() for t in figure.axes[0].texts))
        self.assertFalse(any("q=" in t.get_text() for t in figure.axes[1].texts))

    def test_large_unavailable_run_uses_counts_instead_of_oversized_bitmap(self):
        rows = self.rows(n_loci=60, n_haps=1)
        for row in rows:
            row.update(haplotype_id=None, comparison_status="not_estimable",
                       descriptive_carrier_frequency_difference=None, reason="lead_below_min_call_rate")
        figure = draw_contrasts(self.summary, rows)
        self.assertLess(figure.get_figheight(), 10)
        self.assertTrue(all(not axis.lines and not axis.containers for axis in figure.axes))
        self.assertTrue(any("60 unavailable rows" in text.get_text() for axis in figure.axes for text in axis.texts))

    def test_regional_picture_rejects_multiple_regions(self):
        with self.assertRaisesRegex(ValueError, "exactly one region"):
            draw_region_contrasts(self.summary, self.rows(n_loci=2))


if __name__ == "__main__":
    unittest.main()
