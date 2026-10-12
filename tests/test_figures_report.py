"""Scientific quantities and source labels used by figures and reports."""
from __future__ import annotations

import contextlib
import hashlib
import importlib.util
import io
import json
import math
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parent.parent
CODE = ROOT / "code"
sys.path.insert(0, str(CODE))

import build_report as reports
import compare_haplotypes as pipeline
import make_demo
import make_figures as figures


class FiguresReportTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="test_figures_report_")
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)

    def synthetic_run(self):
        fixtures = make_demo.make_demo(self.folder / "fixtures")
        output = self.folder / "results" / "demo" / "run"
        args = pipeline.argument_parser().parse_args([
            "--vcf", str(fixtures / "cohort.vcf"), "--gwas", str(fixtures / "gwas.tsv"),
            "--samples", str(fixtures / "samples.tsv"), "--out", str(output),
            "--genome-build", "GRCh38", "--gwas-alleles-forward-confirmed",
        ])
        with contextlib.redirect_stdout(io.StringIO()):
            summary = pipeline.run(args)
        return fixtures, output, summary

    def settings(self):
        return {"paths": {
            "out": str(self.folder / "results" / "analysis"),
            "demo_results": str(self.folder / "results" / "demo"),
            "demo_data": str(self.folder / "fixtures"),
        }}

    def test_tiny_p_values_have_finite_correct_negative_log10(self):
        self.assertAlmostEqual(figures.neg_log10_p("3.35e-1003"), 1003 - math.log10(3.35), places=10)
        self.assertAlmostEqual(figures.neg_log10_p("1e-12"), 12.0)
        self.assertEqual(figures.neg_log10_p("1"), 0.0)
        self.assertTrue(math.isfinite(figures.neg_log10_p("3.35e-1003")))
        for invalid in ["0", "-0.1", "1.1", "NaN", "Infinity"]:
            with self.subTest(p=invalid):
                with self.assertRaises(ValueError):
                    figures.neg_log10_p(invalid)

    def test_phase_matrix_uses_resolved_people_over_whole_group(self):
        matrix = figures.prepare_phase_matrix([
            {"locus_id": "L1", "group": "A", "n_resolved": "3", "n_group": "10", "n_called": "9"},
            {"locus_id": "L1", "group": "B", "n_resolved": "0", "n_group": "4", "n_called": "0"},
        ], ["L1"], ["A", "B", "missing"])
        self.assertAlmostEqual(matrix["values"][0][0], 0.3)
        self.assertEqual(matrix["values"][0][1], 0.0)
        self.assertTrue(math.isnan(matrix["values"][0][2]))
        self.assertEqual(matrix["counts"][0][0], "3/10")
        self.assertEqual(matrix["counts"][0][1], "0/4")

    def test_haplotype_matrix_uses_resolved_chromosomes_and_preserves_na(self):
        rows = []
        for group, denominator, copies in [("A", 6, [2, 4]), ("B", 0, [0, 0]), ("C", 6, [0, 6])]:
            for haplotype, count in zip(["H001", "H002"], copies):
                rows.append({
                    "locus_id": "L1", "group": group, "haplotype_id": haplotype,
                    "haplotype_copies": str(count), "resolved_chromosomes": str(denominator),
                    "n_group": "10", "carriers": "1",
                })
        matrix = figures.prepare_haplotype_matrix(rows, "L1", ["A", "B", "C", "missing"])
        h1 = matrix["haplotypes"].index("H001")
        h2 = matrix["haplotypes"].index("H002")
        self.assertAlmostEqual(matrix["values"][h1][0], 2 / 6)
        self.assertAlmostEqual(matrix["values"][h2][0], 4 / 6)
        self.assertTrue(math.isnan(matrix["values"][h1][1]))
        self.assertEqual(matrix["values"][h1][2], 0.0)
        self.assertEqual(matrix["values"][h2][2], 1.0)
        self.assertTrue(math.isnan(matrix["values"][h1][3]))
        self.assertEqual(matrix["denominators"]["A"], 6)

    def test_no_automatic_demo_fallback_and_missing_explicit_run_rejected(self):
        fixtures, _, _ = self.synthetic_run()
        selected, synthetic = reports.select_run(self.settings())
        self.assertIsNone(selected)
        self.assertFalse(synthetic)
        missing = self.folder / "missing_requested_run"
        with self.assertRaises(ValueError):
            reports.select_run(self.settings(), run_dir=missing)
        with self.assertRaises((ValueError, OSError)):
            figures.build_figures(
                gwas_path=fixtures / "gwas.tsv", out_dir=self.folder / "figures",
                run_dir=missing,
            )

    def test_selecting_demo_without_demo_flag_still_identifies_synthetic_data(self):
        _, output, _ = self.synthetic_run()
        selected, synthetic = reports.select_run(self.settings(), run_dir=output)
        self.assertEqual(selected, output)
        self.assertTrue(synthetic)
        configured = self.settings()
        configured["paths"]["out"] = str(output)
        selected, synthetic = reports.select_run(configured)
        self.assertEqual(selected, output)
        self.assertTrue(synthetic)
        configured["paths"]["demo_results"] = str(self.folder / "other_demo_location")
        _, synthetic = reports.select_run(configured, run_dir=output)
        self.assertTrue(synthetic, "Synthetic input provenance must survive an output-directory change")

    def test_run_without_completion_marker_is_rejected_and_skipped_by_demo_selection(self):
        _, completed, _ = self.synthetic_run()
        incomplete = completed.parent / "incomplete_newer_run"
        shutil.copytree(completed, incomplete)
        (incomplete / "COMPLETED.txt").unlink()
        summary_file = incomplete / "run_summary.json"
        summary = json.loads(summary_file.read_text())
        summary["completed_at_utc"] = "2999-01-01T00:00:00+00:00"
        summary_file.write_text(json.dumps(summary))
        with self.assertRaises(ValueError):
            reports.select_run(self.settings(), run_dir=incomplete)
        selected, synthetic = reports.select_run(self.settings(), demo=True)
        self.assertEqual(selected, completed)
        self.assertTrue(synthetic)

    def report_manifest(self, run_dir=None, summary=None):
        return {
            "created_at_utc": "2026-01-01T00:00:00+00:00", "figures": [],
            "sources": {"gwas": str(self.folder / "public_gwas.tsv"), "run_dir": str(run_dir) if run_dir else None},
            "counts": {"gwas_rows": 4, "plotted_rows": 4, "at_or_below_threshold": 3,
                       "above_threshold": 1, "skipped_rows": 0},
            "p_threshold": "1e-8", "run_summary": summary,
        }

    def test_report_keeps_public_gwas_and_unresolved_cohort_status_distinct(self):
        review = {
            "metadata_rows": 4, "vcf_samples": 4, "pasted_category_counts": {"NORM": 2, "T2D": 2},
            "existing_category_counts": {"NORM": 3, "T2D": 1},
            "comparison_counts": {"agreement": 2, "category_conflict": 1, "no_exact_metadata_id": 1},
            "unmatched_metadata_ids": 1, "finalized_sample_table": False,
        }
        content = reports.report_content(self.report_manifest(), review, None, False).lower()
        self.assertIn("public gwas", content)
        self.assertIn("pending", content)
        self.assertIn("metadata", content)
        self.assertIn("conflict", content)

    def test_report_labels_selected_demo_results_as_synthetic(self):
        _, output, summary = self.synthetic_run()
        content = reports.report_content(self.report_manifest(output, summary), None, output, True).lower()
        self.assertIn("synthetic", content)
        self.assertIn("public gwas", content)

    def test_authoritative_groups_are_settled_while_identifiers_are_pending(self):
        review = {
            "authoritative_group_source": True, "finalized_sample_table": False,
            "exact_id_matches": 136, "pending_vcf_id_matches": 11,
            "authoritative_category_counts": {"T2D": 49, "Prediabetes": 47, "NORM": 51},
            "comparison_counts": {"category_conflict": 88},
        }
        manifest = self.report_manifest()
        manifest["figures"] = [{
            "id": "participant_category_counts", "scope": "participant_metadata",
            "title": "Supplied participant categories", "caption": "Source-table counts.",
            "png": str(self.folder / "categories.png"), "pdf": str(self.folder / "categories.pdf"),
        }]
        preflight = {
            "scope": "all_vcf_samples_no_phenotypes", "n_samples": 147, "n_leads": 20,
            "leads_passing_call_rate": 0, "max_high_quality_called": 126,
            "max_high_quality_call_rate": 126 / 147, "raw_reference_genotypes": 0,
            "parameters": {"min_dp": 10, "min_gq": 20, "min_call_rate": 0.9},
        }
        content = reports.report_content(manifest, review, None, False, preflight_summary=preflight)
        self.assertIn("authoritative group source", content)
        self.assertIn("88 historical group-label conflicts", content)
        self.assertIn("resolved in favour of this table", content)
        self.assertIn("group-source choice is settled", content)
        self.assertIn("11 VCF identifier matches await confirmation", content)
        self.assertIn("all 147 VCF samples", content)
        self.assertIn("without using participant group labels", content)
        self.assertIn("No checked lead met", content)
        self.assertIn("126 of 147 (85.7%)", content)
        self.assertLess(content.index("Supplied participant categories"), content.index("Cohort figures are pending"))
        self.assertEqual(content.count("### Supplied participant categories"), 1)
        self.assertNotIn("Final cohort labels remain unresolved", content)
        self.assertNotIn("sample table is finalized", content)
        for group, count in review["authoritative_category_counts"].items():
            self.assertIn(f"| {group} | {count} |", content)

    def test_run_with_no_haplotype_panels_reports_reasons_without_frequency_claims(self):
        output = self.folder / "cohort_run"
        summary = {
            "group_counts": {"T2D": 49, "Prediabetes": 47, "NORM": 51},
            "selected_loci": 20, "locus_status_counts": {"lead_below_min_call_rate": 20},
            "haplotype_carrier_tests_performed": 0,
        }
        content = reports.report_content(self.report_manifest(output, summary), None, output, False)
        self.assertIn("No region yielded an analysable haplotype panel", content)
        self.assertIn("No haplotype frequency estimates were produced", content)
        self.assertIn("Too few participants had a usable lead genotype at the required quality | 20", content)
        self.assertNotIn("Read completeness first, then compare patterns", content)
        self.assertNotIn("Synthetic example", content)

    def test_authoritative_run_claim_requires_matching_saved_sample_table_hash(self):
        fixtures, output, summary = self.synthetic_run()
        supplied = self.folder / "new_authoritative_samples.tsv"
        supplied.write_text("sample\tgroup\nnew_person\tT2D\n")
        review = {
            "authoritative_group_source": True, "finalized_sample_table": True,
            "finalized_sample_table_sha256": hashlib.sha256(supplied.read_bytes()).hexdigest(),
        }
        manifest = self.report_manifest(output, summary)
        content = reports.report_content(manifest, review, output, False)
        self.assertIn("selected run has not been verified", content)
        self.assertNotIn("these group assignments use that table", content)
        # Even another input's matching hash cannot establish sample provenance.
        summary["small_input_sha256"][str(fixtures / "gwas.tsv")] = review["finalized_sample_table_sha256"]
        self.assertIn("selected run has not been verified", reports.report_content(manifest, review, output, False))
        review["finalized_sample_table_sha256"] = hashlib.sha256((fixtures / "samples.tsv").read_bytes()).hexdigest()
        content = reports.report_content(manifest, review, output, False)
        self.assertIn("recorded sample-table checksum matches", content)
        self.assertIn("these group assignments use that table", content)

    def test_render_stages_output_without_exposing_analysis_results_to_quarto_cleanup(self):
        results = self.folder / "results"
        demo_sentinel = results / "demo" / "existing_run" / "COMPLETED.txt"
        figure_sentinel = results / "figures" / "existing_figure.png"
        for path, content in [(demo_sentinel, b"existing completed analysis"),
                              (figure_sentinel, b"existing figure")]:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(content)
        rendered_dirs = []
        expected_html = "<!doctype html><html><body>Staged report</body></html>"

        def fake_quarto(command, *, cwd, check):
            self.assertEqual(Path(cwd), self.folder)
            self.assertTrue(check)
            staging = Path(command[command.index("--output-dir") + 1]).resolve()
            # Real Quarto may clear this entire directory before writing HTML.
            # Fail before writing anything if production sends it our results.
            self.assertNotEqual(staging, results.resolve())
            self.assertFalse(staging.is_relative_to(results.resolve()))
            self.assertTrue(staging.is_relative_to(Path(tempfile.gettempdir()).resolve()))
            self.assertTrue(staging.is_dir())
            rendered_dirs.append(staging)
            (staging / "report.html").write_text(expected_html)

        with mock.patch.object(reports, "PROJECT_ROOT", self.folder), \
                mock.patch.object(reports.shutil, "which", return_value="/fake/quarto"), \
                mock.patch.object(reports.subprocess, "run", side_effect=fake_quarto) as renderer:
            target = reports.render_report()
        renderer.assert_called_once()
        self.assertEqual(target, results / "report.html")
        self.assertEqual(target.read_text(), expected_html)
        self.assertEqual(demo_sentinel.read_bytes(), b"existing completed analysis")
        self.assertEqual(figure_sentinel.read_bytes(), b"existing figure")
        self.assertFalse(rendered_dirs[0].exists())

    @unittest.skipUnless(importlib.util.find_spec("matplotlib"), "Optional matplotlib is not installed")
    def test_optional_export_writes_png_and_pdf_with_synthetic_scope(self):
        fixtures, output, _ = self.synthetic_run()
        manifest = figures.build_figures(
            gwas_path=fixtures / "gwas.tsv", out_dir=self.folder / "figures", run_dir=output,
            synthetic=True, max_loci=1,
        )
        self.assertGreater(len(manifest["figures"]), 0)
        self.assertTrue(any(figure["scope"] == "synthetic_demo" for figure in manifest["figures"]))
        for figure in manifest["figures"]:
            with self.subTest(figure=figure["id"]):
                self.assertTrue(Path(figure["png"]).read_bytes().startswith(b"\x89PNG\r\n\x1a\n"))
                self.assertTrue(Path(figure["pdf"]).read_bytes().startswith(b"%PDF"))


if __name__ == "__main__":
    unittest.main()
