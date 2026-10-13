"""Publication rejects stale evidence and exports one portable aggregate report."""
from __future__ import annotations

from contextlib import redirect_stdout
import csv
import hashlib
from html.parser import HTMLParser
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "code"))
from publish_haplotype_comparison import publish


def table(path, fields, rows=()):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.images, self.links = [], []

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "img":
            self.images.append(attrs["src"])
        elif tag == "a":
            self.links.append(attrs["href"])


class PublicationTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="publication_")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.run, self.evidence, self.out = [self.root / x for x in ("cohort", "evidence", "publication")]
        self.run.mkdir()
        self.evidence.mkdir()
        self.secret = "PRIVATE_PARTICIPANT_84833"
        self.summary = {
            "synthetic": False, "inputs_stable_during_read": True,
            "group_counts": {"T2D": 2, "Prediabetes": 2, "Control": 2},
            "selected_loci": 1, "locus_status_counts": {"lead_below_min_call_rate": 1},
            "haplotype_carrier_tests_performed": 0,
            "parameters": {"test_associations": False, "reference_genotypes_confirmed": False,
                           "min_test_samples": 10, "min_resolved_fraction": .8, "min_carriers": 3,
                           "vcf": str(self.run / (self.secret + ".vcf.gz"))},
            "scan_counts": {"high_quality_homozygous_reference_genotypes": 0},
        }
        self.write_summary()
        (self.run / "COMPLETED.txt").write_text("Completed\n")
        table(self.run / "lead_loci.tsv", ["locus_id", "chrom", "lead_pos", "status"],
              [{"locus_id": "L0001", "chrom": "chr1", "lead_pos": "100", "status": "lead_below_min_call_rate"}])
        table(self.run / "haplotype_group_frequencies.tsv", ["locus_id", "haplotype_id", "group", "frequency"])
        table(self.run / "group_phase_qc.tsv", ["locus_id", "group", "n_resolved"])
        table(self.run / "haplotype_overall_comparisons.tsv", ["locus_id", "test_status", "q_BH_all_overall_tests"])
        table(self.run / "haplotype_comparisons.tsv", ["locus_id", "test_status", "q_BH_all_pairwise_tests"])
        # These private, participant-level artifacts must never join public exports.
        (self.run / "sample_haplotypes.tsv").write_text("sample\thaplotype\n" + self.secret + "\tNA\n")
        (self.run / "samples_used.tsv").write_text("sample\tgroup\n" + self.secret + "\tT2D\n")
        table(self.evidence / "external_gwas_leads.tsv", ["locus_id", "external_gwas_p"],
              [{"locus_id": "L0001", "external_gwas_p": "1e-400"}])
        table(self.evidence / "legacy_locus_audit.tsv", ["legacy_locus", "eligible_for_causal_ranking"])
        table(self.evidence / "candidate_variants.tsv", ["locus_id", "variant", "causal_status"])
        (self.evidence / "REPORT.txt").write_text("No candidates resolved from validated evidence.\n")
        self.audit = {
            "analysis": "aggregate_evidence_audit_not_new_finemapping", "synthetic": False,
            "cohort_group_counts": self.summary["group_counts"], "cohort_selected_loci": 1,
            "cohort_haplotype_association_tests": 0, "resolved_candidate_variants": 0,
            "provenance": {"cohort_summary": {"sha256": digest(self.run / "run_summary.json")},
                           "selected_leads": {"sha256": digest(self.run / "lead_loci.tsv")}},
            "output_sha256": {name: digest(self.evidence / name) for name in
                              ["external_gwas_leads.tsv", "legacy_locus_audit.tsv", "candidate_variants.tsv"]},
        }
        self.write_audit()

    def write_summary(self):
        (self.run / "run_summary.json").write_text(json.dumps(self.summary))

    def write_audit(self):
        (self.evidence / "audit_summary.json").write_text(json.dumps(self.audit))

    def reject_before_export(self):
        with patch("publish_haplotype_comparison.export_panel_c") as exporter:
            with self.assertRaises(ValueError):
                publish(self.run, self.evidence, self.out)
            exporter.assert_not_called()
        self.assertFalse(self.out.exists())

    def test_evidence_from_another_cohort_is_rejected_before_any_output(self):
        self.audit["provenance"]["cohort_summary"]["sha256"] = "0" * 64
        self.write_audit()
        self.reject_before_export()

    def test_modified_lead_table_is_rejected_before_any_output(self):
        path = self.run / "lead_loci.tsv"
        path.write_text(path.read_text().replace("\t100\t", "\t101\t"))
        self.reject_before_export()

    def test_modified_evidence_table_is_rejected_before_any_output(self):
        path = self.evidence / "external_gwas_leads.tsv"
        path.write_text(path.read_text().replace("1e-400", "1e-8"))
        self.reject_before_export()

    def test_synthetic_run_is_rejected_even_with_matching_source_hash(self):
        self.summary["synthetic"] = True
        self.write_summary()
        self.audit["provenance"]["cohort_summary"]["sha256"] = digest(self.run / "run_summary.json")
        self.write_audit()
        self.reject_before_export()

    def test_unstable_run_is_rejected_even_with_matching_source_hash(self):
        self.summary["inputs_stable_during_read"] = False
        self.write_summary()
        self.audit["provenance"]["cohort_summary"]["sha256"] = digest(self.run / "run_summary.json")
        self.write_audit()
        self.reject_before_export()

    def test_incomplete_run_is_rejected_before_any_output(self):
        (self.run / "COMPLETED.txt").unlink()
        self.reject_before_export()

    def test_audit_must_explicitly_identify_real_cohort_evidence(self):
        for synthetic in (True, None):
            with self.subTest(synthetic=synthetic):
                self.audit["synthetic"] = synthetic
                self.write_audit()
                self.reject_before_export()

    def test_evidence_group_counts_must_match_even_when_hashes_are_valid(self):
        self.audit["cohort_group_counts"] = {"T2D": 3, "Prediabetes": 2, "Control": 2}
        self.write_audit()
        self.reject_before_export()

    def test_zero_panels_mean_not_tested_and_export_is_portable_and_aggregate(self):
        def lightweight_export(run_dir, out):
            out.mkdir()
            for name in ["panel_c.png", "panel_c.pdf", "panel_c.svg", "panel_c_data.tsv"]:
                (out / name).write_text("Mocked figure asset\n")
            (out / "REPORT.txt").write_text("Genotype completeness only.\n")
            (out / "report.html").write_text(
                '<html><img src="data:image/png;base64,bW9jay1pbWFnZQ==">'
                '<a href="panel_c.pdf">Figure PDF</a><div class="provenance"></div></html>')
            return {"source_run": str(run_dir.resolve()), "synthetic": False}

        with patch("publish_haplotype_comparison.export_panel_c", side_effect=lightweight_export) as exporter:
            with redirect_stdout(io.StringIO()):
                result = publish(self.run, self.evidence, self.out)
            exporter.assert_called_once()
        self.assertFalse(result["comparison_available"])
        self.assertEqual(result["usable_haplotype_panels"], 0)
        self.assertEqual(result["overall_tests_performed"], 0)
        self.assertEqual(result["pairwise_tests_performed"], 0)
        self.assertEqual(result["candidate_variants_supported_by_cohort_haplotypes"], 0)
        html = (self.out / "report.html").read_text()
        self.assertIn("not tested, not a null result", html)
        self.assertIn("not estimable", html)
        self.assertIn("0 overall and 0 pairwise tests", html)
        self.assertIn("No associated cohort haplotype", html)
        parser = Links()
        parser.feed(html)
        self.assertEqual(len(parser.images), 1)
        self.assertTrue(parser.images[0].startswith("data:image/png;base64,"))
        self.assertNotIn("samples_used.tsv", {p.name for p in self.out.iterdir()})
        self.assertNotIn("sample_haplotypes.tsv", {p.name for p in self.out.iterdir()})
        for path in self.out.iterdir():
            self.assertNotIn(self.secret, path.read_text())
            self.assertNotIn(str(self.root), path.read_text())
        relocated = self.root / "downloaded report"
        self.out.rename(relocated)
        for link in parser.links:
            self.assertFalse(Path(link).is_absolute())
            self.assertTrue((relocated / link).is_file(), link)


if __name__ == "__main__":
    unittest.main()
