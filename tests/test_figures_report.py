"""Scientific quantities and source labels used by figures and reports."""
from __future__ import annotations

import contextlib
import hashlib
import importlib.util
import io
import json
import math
from pathlib import Path
import re
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
from group_labels import display_group, ordered_groups


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

    def test_readable_names_do_not_pool_categories_or_relabel_non_diabetes(self):
        source = {"NORM": 5, "Prediabetes": 3, "T2D": 4, "non_diabetes": 2}
        original = dict(source)
        keys = ordered_groups(source)
        self.assertEqual(keys, ["T2D", "Prediabetes", "NORM", "non_diabetes"])
        self.assertEqual([display_group(key) for key in keys],
                         ["Diabetes", "Prediabetes", "Control", "non_diabetes"])
        self.assertEqual([source[key] for key in keys], [4, 3, 5, 2])
        self.assertEqual(source, original)
        self.assertCountEqual(ordered_groups(["NORM", "Control"]), ["NORM", "Control"])

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

    def test_group_genotype_qc_uses_verified_denominators_and_preserves_missing(self):
        rows = [
            {"locus_id": "L1", "group": "T2D", "n_group": "4", "high_quality_called_samples": "3"},
            {"locus_id": "L1", "group": "Prediabetes", "n_group": "2", "high_quality_called_samples": "0"},
            {"locus_id": "L1", "group": "NORM", "n_group": "0", "high_quality_called_samples": "0"},
            {"locus_id": "L2", "group": "T2D", "n_group": "4", "high_quality_called_samples": "NA"},
        ]
        matrix = figures.prepare_genotype_qc_matrix(rows, ["L1", "L2"], ["T2D", "Prediabetes", "NORM"])
        self.assertEqual(matrix["values"][0][:2], [0.75, 0.0])
        self.assertTrue(math.isnan(matrix["values"][0][2]))
        self.assertTrue(all(math.isnan(value) for value in matrix["values"][1]))
        self.assertEqual(matrix["counts"][0][:2], ["3/4", "0/2"])
        self.assertEqual(matrix["groups"], ["T2D", "Prediabetes", "NORM"])

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
            "exact_id_matches": 58, "pending_vcf_id_matches": 8,
            "authoritative_category_counts": {"T2D": 23, "Prediabetes": 21, "NORM": 22},
            "exact_match_category_counts": {"T2D": 20, "Prediabetes": 19, "NORM": 19},
            "comparison_counts": {"category_conflict": 12},
        }
        manifest = self.report_manifest()
        manifest["figures"] = [{
            "id": "participant_category_counts", "scope": "participant_metadata",
            "title": "Supplied participant categories", "caption": "Source-table counts.",
            "png": str(self.folder / "categories.png"), "pdf": str(self.folder / "categories.pdf"),
        }, {
            "id": "lead_genotype_qc_by_group", "scope": "lead_genotype_qc",
            "title": "Usable lead genotypes by group", "caption": "Exact metadata matches only.",
            "png": str(self.folder / "genotype_qc.png"), "pdf": str(self.folder / "genotype_qc.pdf"),
        }]
        preflight = {
            "scope": "all_vcf_samples_no_phenotypes", "n_samples": 66, "n_leads": 7,
            "leads_passing_call_rate": 0, "max_high_quality_called": 55,
            "max_high_quality_call_rate": 55 / 66, "raw_reference_genotypes": 0,
            "parameters": {"min_dp": 10, "min_gq": 20, "min_call_rate": 0.9},
            "group_qc": {"scope": "exact_metadata_id_matches_only", "matched_sample_count": 58,
                         "unmatched_vcf_count": 8,
                         "matched_category_counts": {"T2D": 20, "Prediabetes": 19, "NORM": 19}},
        }
        content = reports.report_content(manifest, review, None, False, preflight_summary=preflight)
        self.assertIn("authoritative group source", content)
        self.assertIn("12 historical group-label conflicts", content)
        self.assertIn("resolved in favour of this table", content)
        self.assertIn("group-source choice is settled", content)
        self.assertIn("8 VCF identifier matches await confirmation", content)
        self.assertIn("all 66 VCF samples", content)
        self.assertIn("without using participant group labels", content)
        self.assertIn("No checked lead met", content)
        self.assertIn("55 of 66 (83.3%)", content)
        self.assertLess(content.index("Supplied participant categories"), content.index("Cohort haplotype results are unavailable"))
        self.assertEqual(content.count("### Supplied participant categories"), 1)
        self.assertNotIn("Final cohort labels remain unresolved", content)
        self.assertNotIn("sample table is finalized", content)
        self.assertIn("| T2D | Diabetes | 23 | 20 | 3 |", content)
        self.assertIn("| Prediabetes | Prediabetes | 21 | 19 | 2 |", content)
        self.assertIn("| NORM | Control | 22 | 19 | 3 |", content)
        self.assertNotIn("non_diabetes", content)
        self.assertIn("| Diabetes | T2D | 20 |", content)
        self.assertIn("| Prediabetes | Prediabetes | 19 |", content)
        self.assertIn("| Control | NORM | 19 |", content)
        self.assertIn("8 VCF samples remain unassigned", content)
        self.assertIn("does not finalize an analysis subset", content)
        self.assertEqual(content.count("### Usable lead genotypes by group"), 1)

    def test_run_with_no_haplotype_panels_reports_reasons_without_frequency_claims(self):
        output = self.folder / "cohort_run"
        summary = {
            "group_counts": {"T2D": 23, "Prediabetes": 21, "NORM": 22},
            "selected_loci": 7, "locus_status_counts": {"lead_below_min_call_rate": 7},
            "haplotype_carrier_tests_performed": 0,
        }
        content = reports.report_content(self.report_manifest(output, summary), None, output, False)
        self.assertIn("No region yielded an analysable haplotype panel", content)
        self.assertIn("No haplotype frequency estimates were produced", content)
        self.assertIn("Too few participants had a usable lead genotype at the required quality | 7", content)
        self.assertNotIn("Read completeness first, then compare patterns", content)
        self.assertNotIn("Synthetic example", content)

    def test_verified_alias_coverage_is_reported_separately_and_not_left_pending(self):
        counts = {"T2D": 23, "Prediabetes": 21, "Control": 22}
        exact = {"T2D": 20, "Prediabetes": 19, "Control": 19}
        aliases = {"T2D": 3, "Prediabetes": 2, "Control": 3}
        review = {
            "authoritative_group_source": True, "finalized_sample_table": True,
            "metadata_group_column": "Diabetes_Group", "pending_vcf_id_matches": 0,
            "previous_category_change_count": 12, "category_label_only_change_count": 5,
            "authoritative_category_counts": counts, "exact_match_category_counts": exact,
            "verified_alias_category_counts": aliases, "matched_category_counts": counts,
        }
        manifest = self.report_manifest()
        manifest["figures"] = [{
            "id": "lead_genotype_qc_by_group", "scope": "lead_genotype_qc",
            "title": "Usable lead genotypes by group", "caption": "Exact and verified ID matches.",
            "png": str(self.folder / "genotype_qc.png"), "pdf": str(self.folder / "genotype_qc.pdf"),
        }]
        preflight = {
            "n_samples": 66, "n_leads": 7,
            "group_qc": {
                "scope": "exact_and_verified_id_matches", "matched_sample_count": 66,
                "exact_match_sample_count": 58, "verified_alias_sample_count": 8,
                "unmatched_vcf_count": 0, "matched_category_counts": counts,
            },
        }
        content = reports.report_content(manifest, review, None, False, preflight_summary=preflight)
        self.assertIn("`Diabetes_Group` column", content)
        self.assertIn("12 participant group assignments changed", content)
        self.assertIn("5 category labels changed in name only", content)
        self.assertIn("Exact VCF matches | Verified alias matches | Pending ID matches", content)
        self.assertIn("| T2D | Diabetes | 23 | 20 | 3 | 0 |", content)
        self.assertIn("| Prediabetes | Prediabetes | 21 | 19 | 2 | 0 |", content)
        self.assertIn("| Control | Control | 22 | 19 | 3 | 0 |", content)
        self.assertIn("| Total |  | 66 | 58 | 8 | 0 |", content)
        self.assertIn("58 exact metadata–VCF ID matches", content)
        self.assertIn("8 verified alias matches", content)
        self.assertIn("combined exact and verified alias count", content)
        self.assertIn("share of matched participants", content)
        self.assertNotIn("identifier matches still await confirmation", content)
        self.assertNotIn("metadata IDs exactly match the VCF", content)
        self.assertNotIn("share of exact ID matches", content)

        review["matched_category_counts"] = {**counts, "T2D": 22}
        with self.assertRaisesRegex(ValueError, "disagrees with exact and verified alias counts"):
            reports._participant_section(review)

    def test_cohort_overview_is_shown_without_haplotypes_or_saved_preflight(self):
        output = self.folder / "cohort_no_preflight"
        summary = {
            "group_counts": {"T2D": 23, "Prediabetes": 21, "Control": 22},
            "selected_loci": 7,
            "locus_status_counts": {"lead_below_min_call_rate": 6,
                                    "palindromic_lead_orientation_unconfirmed": 1},
            "haplotype_carrier_tests_performed": 0,
            "scan_counts": {"high_quality_homozygous_reference_genotypes": 0},
        }
        manifest = self.report_manifest(output, summary)
        manifest["figures"] = [{
            "id": "cohort_comparison_overview", "scope": "cohort_run",
            "title": "Cohort comparison overview", "caption": "Observed genotypes from saved run tables.",
            "png": str(output / "figures/overview.png"), "pdf": str(output / "figures/overview.pdf"),
        }, {
            "id": "regional_marker_comparison", "scope": "cohort_run",
            "title": "Regional marker comparison", "caption": "Separate SNPs at physical positions.",
            "png": str(output / "figures/regional.png"), "pdf": str(output / "figures/regional.pdf"),
        }]
        content = reports.report_content(manifest, None, output, False,
                                         no_preflight=True, output_dir=output)
        self.assertEqual(content.count("### Cohort comparison overview"), 1)
        self.assertEqual(content.count("### Regional marker comparison"), 1)
        self.assertLess(content.index("No haplotype frequency estimates were produced"),
                        content.index("### Cohort comparison overview"))
        self.assertIn("single-variant summaries, not phased chromosome patterns", content)
        self.assertIn("conditional on observed genotypes", content)
        self.assertIn("no high-quality homozygous-reference (`0/0`) calls", content)
        self.assertIn("ALT allele copies divided by observed chromosome copies", content)
        self.assertIn("No estimate is drawn when its denominator is zero", content)
        self.assertIn("https://www.nature.com/articles/s41588-024-01719-5/figures/1", content)
        self.assertIn("regional variant tables saved by this completed run", content)
        self.assertIn("does not load or rerun a separate VCF precheck", content)
        self.assertIn("No haplotype association tests were run", content)
        self.assertIn("--no-preflight", content)
        self.assertIn("(<figures/overview.png>)", content)
        self.assertNotIn("Lead-genotype completeness check", content)

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

    def test_preflight_and_no_preflight_are_mutually_exclusive(self):
        with mock.patch.object(reports, "load_config") as loader:
            with self.assertRaisesRegex(ValueError, "Choose --preflight or --no-preflight"):
                reports.build_report(preflight=self.folder / "check.json", no_preflight=True)
        loader.assert_not_called()
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as error:
                reports.main(["--preflight", "check.json", "--no-preflight"])
        self.assertEqual(error.exception.code, 2)

    @unittest.skipUnless(importlib.util.find_spec("matplotlib"), "Optional matplotlib is not installed")
    def test_custom_report_skips_saved_preflight_and_preserves_default_and_analysis_files(self):
        fixtures, original_output, _ = self.synthetic_run()
        output = self.folder / "custom report and analysis"
        original_output.rename(output)
        analysis_files = {path: path.read_bytes() for path in output.iterdir() if path.is_file()}
        source = self.folder / "report.qmd"
        source.write_text((ROOT / "report.qmd").read_text())
        config = self.folder / "settings.json"
        metadata = fixtures / "participant_groups.tsv"
        metadata.write_text((fixtures / "samples.tsv").read_text().replace("sample\tgroup\n", "name\tDiabetes_Group\n", 1))
        settings = self.settings()
        settings["paths"].update({"gwas": str(fixtures / "gwas.tsv"),
                                  "metadata": str(metadata)})
        config.write_text(json.dumps(settings))
        defaults = {source: source.read_bytes(), config: config.read_bytes()}
        for name in ["report.html", "report_content.qmd", "report_context.json", "figures/keep.txt"]:
            path = self.folder / "results" / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("Existing default artifact: " + name)
            defaults[path] = path.read_bytes()
        prior_check = self.folder / "data/intermediate/cohort_preflight/summary.json"
        prior_check.parent.mkdir(parents=True)
        prior_check.write_text("Invalid old JSON must never be read when preflight is disabled")
        rendered_dirs = []

        def fake_quarto(command, *, cwd, check):
            self.assertEqual(Path(cwd), output)
            self.assertTrue(check)
            self.assertEqual(command[2], "report.qmd")
            self.assertIn("{{< include report_content.qmd >}}", (output / "report.qmd").read_text())
            content = (output / "report_content.qmd").read_text()
            images = re.findall(r'\]\(<([^>]+)>\)', content)
            downloads = re.findall(r'href="([^"]+\.pdf)"', content)
            self.assertTrue(images)
            self.assertEqual(len(images), len(downloads))
            for relative in images + downloads:
                self.assertTrue(relative.startswith("figures/"), relative)
                self.assertTrue((output / relative).is_file(), relative)
            self.assertNotIn("Lead-genotype completeness check", content)
            self.assertIn("--no-preflight", content)
            self.assertIn("--output-dir", content)
            self.assertIn(str(output), content)
            staging = Path(command[command.index("--output-dir") + 1])
            self.assertFalse(staging.is_relative_to(output))
            rendered_dirs.append(staging)
            shutil.rmtree(staging)
            staging.mkdir()
            (staging / "report.html").write_text("<!doctype html><html>Custom rendered report</html>")

        with mock.patch.object(reports, "PROJECT_ROOT", self.folder), \
                mock.patch.object(reports.shutil, "which", return_value="/fake/quarto"), \
                mock.patch.object(reports.subprocess, "run", side_effect=fake_quarto), \
                mock.patch.object(figures, "build_figures", wraps=figures.build_figures) as builder:
            context = reports.build_report(config=config, run_dir=output, no_preflight=True,
                                           output_dir=output, render=True)
        self.assertIsNone(builder.call_args.kwargs["preflight_path"])
        self.assertIsNone(context["genotype_preflight"])
        self.assertTrue(context["preflight_disabled"])
        self.assertEqual(context["report_source"], str(output / "report.qmd"))
        self.assertEqual(context["html"], str(output / "report.html"))
        stored_context = json.loads((output / "report_context.json").read_text())
        self.assertTrue(stored_context["preflight_disabled"])
        manifest = json.loads(Path(context["figure_manifest"]).read_text())
        self.assertNotIn("lead_genotype_qc", {figure["scope"] for figure in manifest["figures"]})
        for path, expected in {**defaults, **analysis_files}.items():
            self.assertEqual(path.read_bytes(), expected, str(path))
        self.assertFalse(rendered_dirs[0].exists())

    @unittest.skipUnless(importlib.util.find_spec("matplotlib"), "Optional matplotlib is not installed")
    def test_optional_export_writes_png_and_pdf_with_synthetic_scope(self):
        fixtures, output, _ = self.synthetic_run()
        metadata = self.folder / "categories.tsv"
        metadata.write_text(
            "name (Tubingen ID)\tCategory\n"
            "example_a\tNORM\nexample_b\tNORM\nexample_c\tNORM\n"
            "example_d\tPrediabetes\nexample_e\tPrediabetes\nexample_f\tT2D\n"
        )
        qc_table = self.folder / "group_lead_genotype_qc.tsv"
        qc_table.write_text(
            "locus_id\tgroup\tn_group\thigh_quality_called_samples\thigh_quality_call_rate\n"
            "L1\tT2D\t1\t1\t1.0\nL1\tPrediabetes\t2\t1\t0.5\nL1\tNORM\t3\t0\t0.0\n"
            "L2\tT2D\t1\t0\t0.0\nL2\tPrediabetes\t2\t2\t1.0\nL2\tNORM\t3\t3\t1.0\n"
        )
        preflight = self.folder / "preflight.json"
        preflight.write_text(json.dumps({
            "n_samples": 7, "n_leads": 2,
            "parameters": {"min_dp": 10, "min_gq": 20, "min_call_rate": 0.9},
            "inputs": {"metadata": {"path": str(metadata), "sha256": hashlib.sha256(metadata.read_bytes()).hexdigest()}},
            "group_qc": {
                "scope": "exact_metadata_id_matches_only", "table": str(qc_table),
                "source_category_counts": {"T2D": 1, "Prediabetes": 2, "NORM": 3},
                "matched_category_counts": {"T2D": 1, "Prediabetes": 2, "NORM": 3},
                "matched_sample_count": 6, "unmatched_vcf_count": 1, "unmatched_metadata_count": 0,
            },
        }))
        with mock.patch.object(figures, "_save", wraps=figures._save) as saved:
            manifest = figures.build_figures(
                gwas_path=fixtures / "gwas.tsv", out_dir=self.folder / "figures", run_dir=output,
                synthetic=True, max_loci=1, metadata_path=metadata, preflight_path=preflight,
            )
        self.assertGreater(len(manifest["figures"]), 0)
        self.assertTrue(any(figure["scope"] == "synthetic_demo" for figure in manifest["figures"]))
        self.assertEqual(manifest["group_label_map"],
                         {"T2D": "Diabetes", "Prediabetes": "Prediabetes", "NORM": "Control"})
        self.assertEqual(manifest["counts"]["participant_category_counts"],
                         {"NORM": 3, "Prediabetes": 2, "T2D": 1})
        plot_calls = {call.args[3]: call.args[1] for call in saved.call_args_list}
        group_labels = ["Diabetes", "Prediabetes", "Control"]
        for identifier, plot in plot_calls.items():
            if (identifier in {"participant_category_counts", "lead_genotype_qc_by_group", "phase_resolution"}
                    or identifier.startswith("haplotype_frequencies_")):
                self.assertEqual([tick.get_text() for tick in plot.axes[0].get_xticklabels()], group_labels)
        category_axes = plot_calls["participant_category_counts"].axes[0]
        self.assertEqual([bar.get_height() for bar in category_axes.patches], [1, 2, 3])
        qc_axes = plot_calls["lead_genotype_qc_by_group"].axes[0]
        self.assertEqual(qc_axes.images[0].get_array().tolist(), [[1.0, 0.5, 0.0], [0.0, 1.0, 1.0]])
        qc_figure = next(item for item in manifest["figures"] if item["scope"] == "lead_genotype_qc")
        self.assertNotIn("SYNTHETIC DEMO", qc_figure["title"])
        # The new layouts place group labels on the count panel and regional rows.
        # They must describe the synthetic run's four people per group, rather
        # than reuse the separately supplied metadata's source-table counts.
        overview = plot_calls["cohort_comparison_overview"]
        self.assertEqual([tick.get_text() for tick in overview.axes[1].get_xticklabels()], group_labels)
        self.assertEqual([bar.get_height() for bar in overview.axes[1].patches], [4, 4, 4])
        self.assertEqual([label.get_text() for label in overview.axes[2].get_legend().get_texts()],
                         [*group_labels, "All participants"])
        regional = plot_calls["regional_marker_comparison"]
        self.assertEqual([axis.get_ylabel() for axis in regional.axes[::2]],
                         [f"{group}\nn = 4\nPercent" for group in group_labels])
        for identifier in ["cohort_comparison_overview", "regional_marker_comparison"]:
            item = next(figure for figure in manifest["figures"] if figure["id"] == identifier)
            self.assertEqual(item["scope"], "synthetic_demo")
            self.assertTrue(item["title"].startswith("SYNTHETIC DEMO"))
            self.assertTrue(item["caption"].startswith("SYNTHETIC DEMO"))
            self.assertTrue(any(text.get_text().startswith("SYNTHETIC DEMO")
                                for text in plot_calls[identifier].texts))
        for figure in manifest["figures"]:
            with self.subTest(figure=figure["id"]):
                self.assertTrue(Path(figure["png"]).read_bytes().startswith(b"\x89PNG\r\n\x1a\n"))
                self.assertTrue(Path(figure["pdf"]).read_bytes().startswith(b"%PDF"))
        # A different category source must not silently reuse old group QC.
        metadata.write_text(metadata.read_text().replace("example_f\tT2D", "example_f\tNORM"))
        with self.assertRaisesRegex(ValueError, "(?i)(metadata|category|checksum)"):
            figures.build_figures(
                gwas_path=fixtures / "gwas.tsv", out_dir=self.folder / "stale_figures",
                metadata_path=metadata, preflight_path=preflight,
            )


if __name__ == "__main__":
    unittest.main()
