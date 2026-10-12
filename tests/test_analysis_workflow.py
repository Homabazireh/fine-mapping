"""Output and completion invariants for the five-step analysis workflow."""
from __future__ import annotations

import contextlib
import csv
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "code"))

import compare_haplotypes as pipeline
from make_demo import make_demo


class AnalysisWorkflowTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="test_analysis_workflow_")
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)
        self.fixtures = make_demo(self.folder / "inputs")
        self.out = self.folder / "results"
        self.args = pipeline.argument_parser().parse_args([
            "--vcf", str(self.fixtures / "cohort.vcf"),
            "--gwas", str(self.fixtures / "gwas.tsv"),
            "--samples", str(self.fixtures / "samples.tsv"),
            "--out", str(self.out), "--genome-build", "GRCh38",
        ])

    def read_output(self, filename):
        with (self.out / filename).open() as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def run_analysis(self):
        with contextlib.redirect_stdout(io.StringIO()):
            return pipeline.run(self.args)

    def test_group_tables_conserve_resolved_people_and_chromosome_copies(self):
        summary = self.run_analysis()
        samples = self.read_output("sample_haplotypes.tsv")
        frequencies = self.read_output("haplotype_group_frequencies.tsv")
        group_qc = self.read_output("group_phase_qc.tsv")
        self.assertEqual(summary["group_counts"], {"T2D": 4, "Prediabetes": 4, "NORM": 4})
        self.assertEqual(summary["locus_status_counts"], {"analysed": 2})
        self.assertEqual(summary["haplotype_carrier_tests_performed"], 0)
        self.assertEqual(len(self.read_output("selected_markers.tsv")), 7)

        assignments = {(row["locus_id"], row["sample"]): row for row in samples}
        self.assertEqual(assignments["L0001", "D4"]["status"], "different_phase_sets")
        self.assertEqual(assignments["L0001", "N4"]["status"], "missing_or_low_quality")
        for sample in ["D4", "N4"]:
            self.assertEqual(assignments["L0001", sample]["haplotype_1"], "NA")
            self.assertEqual(assignments["L0001", sample]["haplotype_2"], "NA")

        for qc in group_qc:
            with self.subTest(locus=qc["locus_id"], group=qc["group"]):
                group_samples = [row for row in samples
                                 if row["locus_id"] == qc["locus_id"] and row["group"] == qc["group"]]
                resolved = [row for row in group_samples if row["haplotype_1"] != "NA"]
                rows = [row for row in frequencies
                        if row["locus_id"] == qc["locus_id"] and row["group"] == qc["group"]]
                self.assertEqual(int(qc["n_group"]), len(group_samples))
                self.assertEqual(int(qc["n_resolved"]), len(resolved))
                self.assertEqual(int(qc["n_unresolved"]), len(group_samples) - len(resolved))
                self.assertEqual(sum(int(row["haplotype_copies"]) for row in rows), 2 * len(resolved))
                for row in rows:
                    haplotype = row["haplotype_id"]
                    pairs = [[sample["haplotype_1"], sample["haplotype_2"]] for sample in resolved]
                    copies = sum(pair.count(haplotype) for pair in pairs)
                    carriers = sum(haplotype in pair for pair in pairs)
                    self.assertEqual(int(row["haplotype_copies"]), copies)
                    self.assertEqual(int(row["carriers"]), carriers)
                    self.assertEqual(int(row["resolved_chromosomes"]), 2 * len(resolved))
                    self.assertAlmostEqual(float(row["haplotype_frequency"]), copies / (2 * len(resolved)), places=7)
                    self.assertAlmostEqual(float(row["carrier_frequency"]), carriers / len(resolved), places=7)

        self.assertEqual((self.out / "regions.bed").read_text().splitlines(),
                         ["chr1\t0\t60000", "chr1\t150009\t250010"])
        self.assertTrue(summary["inputs_stable_during_read"])
        self.assertTrue((self.out / "COMPLETED.txt").is_file())

    def test_changed_input_during_scan_cannot_receive_completion_marker(self):
        original_scan = pipeline.scan_vcf

        def scan_then_change_input(*args, **kwargs):
            counts = original_scan(*args, **kwargs)
            with (self.fixtures / "cohort.vcf").open("a") as handle:
                handle.write("# changed after the VCF scan\n")
            return counts

        with mock.patch.object(pipeline, "scan_vcf", side_effect=scan_then_change_input):
            with self.assertRaisesRegex(ValueError, "input changed during the run"):
                self.run_analysis()
        self.assertFalse((self.out / "COMPLETED.txt").exists())
        summary = json.loads((self.out / "run_summary.json").read_text())
        self.assertFalse(summary["inputs_stable_during_read"])

    def test_unresolved_sample_groups_fail_before_creating_outputs(self):
        samples = self.fixtures / "samples.tsv"
        samples.write_text(samples.read_text().replace("D1\tT2D", "D1\tUNRESOLVED"))
        with self.assertRaisesRegex(ValueError, "Unresolved sample/group"):
            self.run_analysis()
        self.assertFalse(self.out.exists())

    def test_existing_analysis_outputs_cannot_be_overwritten(self):
        self.out.mkdir()
        sentinel = self.out / "existing_results.txt"
        sentinel.write_text("keep previous analysis")
        with self.assertRaisesRegex(ValueError, "Output directory is not empty"):
            self.run_analysis()
        self.assertEqual(sentinel.read_text(), "keep previous analysis")
        self.assertEqual(list(self.out.iterdir()), [sentinel])


if __name__ == "__main__":
    unittest.main()
