"""The private contrast figure must not turn missing comparisons into zero effects."""
from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "code"))

from compare_haplotypes import COMP_FIELDS, FREQ_FIELDS
from plot_haplotype_contrasts import INFERENCE_SOURCE_FIELDS, draw_contrasts, export_contrasts, interpretation, load_contrasts


def table(path, fields, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "NA") for key in fields})


class HaplotypeContrastTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="contrast_figure_test_")
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)
        self.run = self.folder / "run"
        self.run.mkdir()
        self.summary = {"synthetic": True, "inputs_stable_during_read": True,
                        "group_counts": {"T2D": 24, "Prediabetes": 23, "Control": 22},
                        "selected_loci": 2,
                        "locus_status_counts": {"lead_below_min_call_rate": 1, "palindromic_lead_orientation_unconfirmed": 1},
                        "completed_at_utc": "2026-10-13T00:00:00+00:00",
                        "parameters": {"test_associations": True, "reference_genotypes_confirmed": False,
                                       "min_call_rate": .9, "min_test_samples": 10, "min_resolved_fraction": .8,
                                       "min_carriers": 3}}
        self.leads = [{"locus_id": "L0001", "chrom": "chr1", "lead_pos": "100", "status": "lead_below_min_call_rate"},
                      {"locus_id": "L0002", "chrom": "chr2", "lead_pos": "200", "status": "palindromic_lead_orientation_unconfirmed"}]
        self.save()

    def save(self, pairs=(), frequencies=()):
        (self.run / "COMPLETED.txt").write_text("Completed synthetic fixture")
        (self.run / "run_summary.json").write_text(json.dumps(self.summary))
        table(self.run / "lead_loci.tsv", list(self.leads[0]), self.leads)
        table(self.run / "haplotype_comparisons.tsv", COMP_FIELDS, pairs)
        table(self.run / "haplotype_group_frequencies.tsv", FREQ_FIELDS, frequencies)

    def eligible_fixture(self):
        self.leads[0]["status"] = "analysed"
        self.summary["locus_status_counts"] = {"analysed": 1, "palindromic_lead_orientation_unconfirmed": 1}
        self.summary["parameters"]["reference_genotypes_confirmed"] = True
        pair = {"locus_id": "L0001", "haplotype_id": "L0001_H001", "group_a": "Control", "group_b": "T2D",
                "n_resolved_a": 20, "n_resolved_b": 20, "carriers_a": 6, "carriers_b": 10,
                "carrier_frequency_a": .3, "carrier_frequency_b": .5, "test_status": "tested_unadjusted",
                "carrier_frequency_difference_a_minus_b": -.2, "carrier_difference_CI95_low": -.4,
                "carrier_difference_CI95_high": -.02, "p_fisher_two_sided": .01,
                "q_BH_all_pairwise_tests": .03, "q_BH_all_haplotype_tests": .03,
                "conditional_odds_ratio_a_vs_b": .428571, "CI95_low": .1, "CI95_high": .9}
        frequencies = [{"locus_id": "L0001", "haplotype_id": "L0001_H001", "group": group, "haplotype_frequency": value}
                       for group, value in [("T2D", .3), ("Prediabetes", .2), ("Control", .18)]]
        self.save([pair], frequencies)
        return pair, frequencies

    def test_failed_panels_have_both_contrasts_and_no_invented_haplotypes(self):
        summary, rows = load_contrasts(self.run)
        self.assertEqual(len(rows), 4)
        self.assertEqual({r["contrast"] for r in rows}, {"Diabetes vs Control", "Prediabetes vs Control"})
        for row in rows:
            self.assertEqual(row["comparison_status"], "not_estimable")
            self.assertEqual(row["reason"], row["panel_status"])
            for field in ["haplotype_id", "carrier_frequency_difference", "CI95_low", "CI95_high", "q_BH_all_pairwise_tests"]:
                self.assertIsNone(row[field])
        text = interpretation(summary, rows)
        self.assertIn("0 requested haplotype comparisons were tested", text)
        self.assertIn("does not mean zero effect", text)
        self.assertIn("No candidate causal variant can be localized", text)

    def test_reversed_pair_flips_difference_and_interval_without_recorrecting_q(self):
        self.eligible_fixture()
        _, rows = load_contrasts(self.run)
        row = next(r for r in rows if r["comparison_status"] == "tested")
        self.assertEqual(row["contrast"], "Diabetes vs Control")
        self.assertEqual((row["carrier_frequency_difference"], row["CI95_low"], row["CI95_high"]), (.2, .02, .4))
        self.assertEqual((row["carriers_case"], row["carriers_control"]), (10, 6))
        self.assertEqual((row["n_resolved_case"], row["n_resolved_control"]), (20, 20))
        self.assertEqual(row["q_BH_all_pairwise_tests"], .03)
        self.assertEqual(row["haplotype_frequency_case_descriptive"], .3)
        self.assertNotEqual(row["haplotype_frequency_case_descriptive"], row["carrier_frequency_case"])
        pred = next(r for r in rows if r["locus_id"] == "L0001" and r["contrast"] == "Prediabetes vs Control")
        self.assertEqual(pred["reason"], "pairwise_result_unavailable")

    def test_skipped_haplotype_does_not_plot_a_descriptive_difference_as_inference(self):
        pair, freq = self.eligible_fixture()
        pair["test_status"] = "not_tested_reference_provenance"
        for field in INFERENCE_SOURCE_FIELDS:
            pair[field] = "NA"
        self.summary["parameters"]["reference_genotypes_confirmed"] = False
        self.save([pair], freq)
        _, rows = load_contrasts(self.run)
        row = next(r for r in rows if r["reason"] == "not_tested_reference_provenance")
        self.assertIsNone(row["carrier_frequency_difference"])
        self.assertIsNone(row["q_BH_all_pairwise_tests"])
        self.assertEqual(row["haplotype_frequency_case_descriptive"], .3)
        self.assertEqual(row["comparison_status"], "descriptive_only")
        self.assertEqual(row["inference_status"], "not_tested_reference_provenance")
        self.assertAlmostEqual(row["descriptive_carrier_frequency_difference"], .2)

    def descriptive_fixture(self):
        pair, freq = self.eligible_fixture()
        pair["test_status"] = "not_tested_reference_provenance"
        for field in INFERENCE_SOURCE_FIELDS:
            pair[field] = "NA"
        self.summary["parameters"]["reference_genotypes_confirmed"] = False
        return pair, freq

    def test_true_observed_zero_is_distinct_from_missing_and_has_no_inference(self):
        pair, freq = self.descriptive_fixture()
        pair.update(carriers_a=10, carrier_frequency_a=.5)
        self.save([pair], freq)
        _, rows = load_contrasts(self.run)
        observed = rows[0]
        self.assertEqual(observed["descriptive_carrier_frequency_difference"], 0)
        self.assertEqual(observed["comparison_status"], "descriptive_only")
        self.assertIsNone(observed["CI95_low"])
        self.assertIsNone(rows[1]["descriptive_carrier_frequency_difference"])
        self.assertEqual(rows[1]["comparison_status"], "not_estimable")

    def test_untested_counts_frequencies_and_inference_fields_are_validated(self):
        pair, freq = self.descriptive_fixture()
        for update, expected in [({"carriers_b": 21}, "counts"),
                                 ({"carriers_b": 9.5}, "counts"),
                                 ({"carrier_frequency_b": .6}, "disagrees"),
                                 ({"p_fisher_two_sided": .01}, "Untested")]:
            with self.subTest(update=update):
                self.save([dict(pair, **update)], freq)
                with self.assertRaisesRegex(ValueError, expected):
                    load_contrasts(self.run)

    def test_zero_denominator_stays_unavailable(self):
        pair, freq = self.descriptive_fixture()
        pair.update(n_resolved_b=0, carriers_b=0, carrier_frequency_b="NA")
        self.save([pair], freq)
        _, rows = load_contrasts(self.run)
        self.assertIsNone(rows[0]["descriptive_carrier_frequency_difference"])
        self.assertEqual(rows[0]["comparison_status"], "not_estimable")

    def test_descriptive_difference_uses_exact_counts_not_rounded_frequency(self):
        pair, freq = self.descriptive_fixture()
        pair.update(n_resolved_b=21, carriers_b=7, carrier_frequency_b=.33333333)
        self.save([pair], freq)
        _, rows = load_contrasts(self.run)
        self.assertEqual(rows[0]["descriptive_carrier_frequency_difference"], 7/21 - 6/20)

    def test_tested_rows_cannot_bypass_provenance_or_resolution_gates(self):
        pair, freq = self.eligible_fixture()
        self.summary["parameters"]["reference_genotypes_confirmed"] = False
        self.save([pair], freq)
        with self.assertRaisesRegex(ValueError, "reference-call provenance"):
            load_contrasts(self.run)
        self.summary["parameters"]["reference_genotypes_confirmed"] = True
        pair["n_resolved_b"] = 10
        pair["carriers_b"] = 5
        self.save([pair], freq)
        with self.assertRaisesRegex(ValueError, "resolution"):
            load_contrasts(self.run)

    def test_duplicate_and_contradictory_effect_rows_fail(self):
        pair, freq = self.eligible_fixture()
        self.save([pair, pair], freq)
        with self.assertRaisesRegex(ValueError, "Duplicate pairwise"):
            load_contrasts(self.run)
        pair["carrier_frequency_difference_a_minus_b"] = -.1
        self.save([pair], freq)
        with self.assertRaisesRegex(ValueError, "disagrees"):
            load_contrasts(self.run)

    def test_incomplete_or_unstable_run_cannot_be_rendered(self):
        (self.run / "COMPLETED.txt").unlink()
        with self.assertRaisesRegex(ValueError, "completed"):
            load_contrasts(self.run)
        self.summary["inputs_stable_during_read"] = False
        self.save()
        with self.assertRaisesRegex(ValueError, "stability"):
            load_contrasts(self.run)

    @unittest.skipUnless(importlib.util.find_spec("matplotlib"), "Matplotlib optional")
    def test_unavailable_figure_has_no_effect_points_or_zero_lines(self):
        from matplotlib import pyplot as plt
        summary, rows = load_contrasts(self.run)
        figure = draw_contrasts(summary, rows)
        self.addCleanup(plt.close, figure)
        for axis in figure.axes:
            self.assertEqual(len(axis.lines), 0)
            self.assertEqual(len(axis.collections), 0)
            self.assertEqual(len(axis.get_xticks()), 0)
            self.assertTrue(all("Not estimable" in t.get_text() for t in axis.texts))

    @unittest.skipUnless(importlib.util.find_spec("matplotlib"), "Matplotlib optional")
    def test_export_is_private_and_preserves_existing_artifacts(self):
        out = self.folder / "contrasts"
        manifest = export_contrasts(self.run, out)
        self.assertTrue(manifest["private_results"])
        self.assertEqual(manifest["comparisons_tested"], 0)
        self.assertEqual(manifest["comparisons_not_estimable"], 4)
        for suffix in ["png", "pdf", "svg"]:
            self.assertTrue((out / f"haplotype_contrasts.{suffix}").is_file())
        with (out / "contrast_data.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertTrue(all(r["carrier_frequency_difference"] == "NA" for r in rows))
        with self.assertRaisesRegex(ValueError, "new or empty"):
            export_contrasts(self.run, out)


if __name__ == "__main__":
    unittest.main()
