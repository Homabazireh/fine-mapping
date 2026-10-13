"""Descriptive copy-frequency displays preserve complete local patterns and missingness."""
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
from plot_descriptive_haplotypes import export_descriptive, load_descriptive, region_interpretations


class DescriptiveHaplotypeTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="descriptive_haplotype_fixture_")
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)
        self.run = self.folder / "run"
        self.run.mkdir()
        (self.run / "COMPLETED.txt").write_text("Artificial completed fixture")
        self.summary = {"synthetic": True, "inputs_stable_during_read": True,
                        "selected_loci": 2, "distance_clusters_available": 2,
                        "group_counts": {"T2D": 5, "Prediabetes": 6, "Control": 4},
                        "haplotype_carrier_tests_performed": 0, "haplotype_overall_tests_performed": 0}
        self.leads = [{"locus_id": "L0001", "chrom": "chr1", "lead_pos": "100", "status": "analysed"},
                      {"locus_id": "L0002", "chrom": "chr2", "lead_pos": "200", "status": "lead_below_min_call_rate"}]
        self.rows = []
        for group, resolved, counts in [("T2D", 5, [5, 3, 2]), ("Prediabetes", 0, [0, 0, 0]), ("Control", 4, [2, 4, 2])]:
            for index, copies in enumerate(counts, 1):
                self.rows.append({"locus_id": "L0001", "haplotype_id": f"L0001_H{index:03d}", "group": group,
                                  "n_group": self.summary["group_counts"][group], "n_resolved": resolved,
                                  "haplotype_copies": copies, "resolved_chromosomes": 2*resolved,
                                  "haplotype_frequency": f"{copies/(2*resolved):.8g}" if resolved else "NA"})
        self.save()

    def save(self):
        (self.run / "run_summary.json").write_text(json.dumps(self.summary))
        for name, rows in [("lead_loci.tsv", self.leads), ("haplotype_group_frequencies.tsv", self.rows)]:
            with (self.run / name).open("w", newline="") as handle:
                writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
                writer.writeheader(); writer.writerows(rows)

    def test_all_patterns_use_copy_counts_and_both_control_contrasts(self):
        _, plotted, differences = load_descriptive(self.run)
        self.assertEqual(len(plotted), 9)
        self.assertEqual(len(differences), 6)
        first = next(r for r in differences if r["haplotype_id"] == "L0001_H001" and r["contrast"] == "Diabetes vs Control")
        self.assertEqual(first["haplotype_frequency_case"], .5)
        self.assertEqual(first["haplotype_frequency_control"], .25)
        self.assertEqual(first["descriptive_copy_frequency_difference"], .25)
        self.assertEqual(first["resolved_chromosomes_case"], 10)
        self.assertEqual(first["n_resolved_case"], 5)
        self.assertTrue(all(r["interpretation_status"] == "descriptive_only_not_tested" for r in differences))
        self.assertFalse(any("p_value" in r or "q_value" in r or "CI95_low" in r for r in differences))

    def test_missing_group_stays_unavailable_not_zero_frequency(self):
        _, plotted, differences = load_descriptive(self.run)
        unavailable = [r for r in plotted if r["group"] == "Prediabetes"]
        self.assertTrue(all(r["haplotype_frequency"] is None and r["resolved_chromosomes"] == 0 for r in unavailable))
        self.assertTrue(all(r["descriptive_copy_frequency_difference"] is None for r in differences if r["case_group"] == "Prediabetes"))
        notes = region_interpretations(differences)
        pred = next(r for r in notes if r["contrast"] == "Prediabetes vs Control")
        self.assertIsNone(pred["haplotype_id_largest_absolute_descriptive_difference"])

    def test_largest_descriptive_ties_are_stable_and_never_significance_claims(self):
        _, _, differences = load_descriptive(self.run)
        selected = [r for r in differences if r["contrast"] == "Diabetes vs Control"]
        selected[0]["descriptive_copy_frequency_difference"] = .2
        selected[1]["descriptive_copy_frequency_difference"] = -.20000000000000004
        note = region_interpretations(selected)[0]
        self.assertEqual(note["haplotype_id_largest_absolute_descriptive_difference"], "L0001_H001")
        self.assertIn("Descriptive only", note["interpretation"])
        self.assertIn("No association test", note["interpretation"])

    def test_bad_copy_totals_and_missing_group_rows_are_rejected(self):
        self.rows[0]["haplotype_copies"] = 4
        self.rows[0]["haplotype_frequency"] = .4
        self.save()
        with self.assertRaisesRegex(ValueError, "do not sum"):
            load_descriptive(self.run)
        self.rows.pop(0)
        self.save()
        with self.assertRaisesRegex(ValueError, "all three group rows"):
            load_descriptive(self.run)

    def test_scope_completed_run_and_no_inference_requirements(self):
        self.summary["distance_clusters_available"] = 3
        self.save()
        with self.assertRaisesRegex(ValueError, "all available"):
            load_descriptive(self.run)
        self.summary["distance_clusters_available"] = 2
        self.summary["haplotype_overall_tests_performed"] = 1
        self.save()
        with self.assertRaisesRegex(ValueError, "no inferential tests"):
            load_descriptive(self.run)
        (self.run / "COMPLETED.txt").unlink()
        with self.assertRaisesRegex(ValueError, "completed"):
            load_descriptive(self.run)

    @unittest.skipUnless(importlib.util.find_spec("matplotlib"), "Matplotlib optional")
    def test_exports_single_figure_all_formats_and_private_tables(self):
        out = self.folder / "figures"
        manifest = export_descriptive(self.run, out, facets_per_row=2, dpi=72)
        self.assertEqual(manifest["local_patterns"], 3)
        self.assertEqual(manifest["contrast_rows"], 6)
        self.assertTrue(manifest["private_results_only"])
        self.assertTrue(manifest["descriptive_only"])
        for extension in ["png", "pdf", "svg"]:
            self.assertTrue((out / f"descriptive_haplotypes.{extension}").is_file())
        self.assertIn("No inferential association tests were performed", (out / "INTERPRETATION.txt").read_text())
        with self.assertRaisesRegex(ValueError, "new or empty"):
            export_descriptive(self.run, out)


if __name__ == "__main__":
    unittest.main()
