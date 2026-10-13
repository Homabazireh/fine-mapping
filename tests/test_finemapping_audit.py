"""Fine-mapping audits must preserve evidence boundaries and original inputs."""
from __future__ import annotations

import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "code"))
from audit_finemapping import (
    alignment_errors, audit_finemapping, external_leads, retained_genotype_counts,
)


def table(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class FineMappingAuditTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.run, self.legacy, self.out = [self.root / name for name in ("cohort", "legacy", "audit")]
        self.run.mkdir()
        (self.run / "COMPLETED.txt").write_text("Complete\n")
        self.summary = {
            "genome_build_user_declared": "GRCh38", "group_counts": {"Control": 2, "T2D": 2, "Prediabetes": 2},
            "haplotype_carrier_tests_performed": 0,
            "scan_counts": {"high_quality_homozygous_reference_genotypes": 0},
        }
        (self.run / "run_summary.json").write_text(json.dumps(self.summary))
        self.lead = {"locus_id": "L0001", "chrom": "chr1", "lead_pos": "12", "lead_effect_allele": "A",
                     "lead_other_allele": "G", "lead_p": "3.35E-1003", "status": "lead_below_min_call_rate",
                     "selected_markers": "0"}
        table(self.run / "lead_loci.tsv", list(self.lead), [self.lead])
        self.gwas = self.root / "gwas.tsv"
        self.gwas_row = {"Chromsome": "1", "Position": "12", "EffectAllele": "A", "NonEffectAllele": "G",
                         "Pval": "3.35e-1003", "Beta": "-0.2", "SE": "0.01", "Neff": "12345",
                         "liftover_strand": "+"}
        table(self.gwas, list(self.gwas_row), [self.gwas_row])

    def create_legacy(self):
        table(self.legacy / "all_loci_summary.tsv", ["CHR", "LEAD_BP", "STATUS"],
              [{"CHR": "chr1", "LEAD_BP": "12", "STATUS": "SUCCESS"}])
        locus = self.legacy / "loci" / "chr1_12"
        qc = {"ONT_samples": "6", "Variants_analyzed": "1", "Converged": "TRUE", "GWAS_n": "NOT_SUPPLIED",
              "LD_shrinkage": "0.02"}
        table(locus / "susie_QC.tsv", list(qc), [qc])
        pip = {"Variant_key": "chr1:12:A:G", "GWAS_EA": "A", "GWAS_NEA": "G", "ONT_REF": "A", "ONT_ALT": "G",
               "Beta": "-0.2", "SE": "0.01", "Z_ALT": "20", "PIP": "0.99999123456789"}
        table(locus / "susie_PIP.tsv", list(pip), [pip])
        (locus / "locus_genotypes.tsv").write_text("chr1\t12\tA\tG\t0|1\t1/0\t1|1\t./.\t.|.\t1/1\n")
        return locus, pip

    def test_tiny_pvalues_effect_sign_and_missing_sample_sizes_are_preserved(self):
        row = external_leads([self.lead], self.gwas)[0]
        self.assertEqual(row["external_gwas_p"], "3.35e-1003")
        self.assertEqual(row["external_gwas_neff"], "12345")
        self.assertEqual(row["external_gwas_ncases"], "NA")
        self.assertEqual(row["effect_direction"], "negative")
        self.assertEqual(row["causal_status"], "unresolved_external_association_only")

    def test_high_pip_and_convergence_do_not_become_causal_claims(self):
        locus, pip = self.create_legacy()
        original = (locus / "locus_genotypes.tsv").read_bytes()
        result = audit_finemapping(self.run, self.gwas, self.legacy, self.out)
        self.assertEqual(result["resolved_candidate_variants"], 0)
        self.assertEqual(result["legacy_eligible_for_causal_ranking"], 0)
        self.assertEqual(result["legacy_exclusion_counts"]["no_reference_homozygotes_in_retained_ld_genotypes"], 1)
        self.assertEqual(result["legacy_exclusion_counts"]["gwas_sample_size_not_supplied"], 1)
        self.assertEqual((locus / "locus_genotypes.tsv").read_bytes(), original)
        for path in self.out.iterdir():
            self.assertNotIn(pip["PIP"], path.read_text())
            self.assertNotIn(str(self.root), path.read_text())
        self.assertEqual(len((self.out / "candidate_variants.tsv").read_text().splitlines()), 1)
        with self.assertRaisesRegex(ValueError, "new or empty"):
            audit_finemapping(self.run, self.gwas, self.legacy, self.out)

    def test_missing_genotypes_are_never_treated_as_reference(self):
        locus, _ = self.create_legacy()
        counts = retained_genotype_counts(locus / "locus_genotypes.tsv", {"chr1:12:A:G"}, 6)
        self.assertEqual(counts["ref_hom"], 0)
        self.assertEqual(counts["het"], 2)
        self.assertEqual(counts["alt_hom"], 2)
        self.assertEqual(counts["missing_or_other"], 2)
        self.assertEqual(counts["width_mismatch_rows"], 0)
        self.assertEqual(retained_genotype_counts(locus / "locus_genotypes.tsv", {"chr1:12:A:G"}, 7)["width_mismatch_rows"], 1)

    def test_alignment_check_detects_ref_flip_error(self):
        _, row = self.create_legacy()
        self.assertEqual(alignment_errors([row]), (0, 0))
        self.assertEqual(alignment_errors([{**row, "Z_ALT": "-20"}]), (1, 0))
        self.assertEqual(alignment_errors([{**row, "GWAS_NEA": "C"}]), (1, 0))

    def test_source_mismatch_and_duplicate_leads_fail_before_outputs(self):
        for rows in ([{**self.gwas_row, "Pval": "1e-9"}], [self.gwas_row, self.gwas_row],
                     [{**self.gwas_row, "EffectAllele": "G", "NonEffectAllele": "A"}]):
            table(self.gwas, list(self.gwas_row), rows)
            with self.assertRaises(ValueError):
                audit_finemapping(self.run, self.gwas, self.legacy, self.out)
            self.assertFalse(self.out.exists())

    def test_no_legacy_artifact_does_not_invent_evidence(self):
        result = audit_finemapping(self.run, self.gwas, self.legacy, self.out)
        self.assertEqual(result["legacy_loci"], 0)
        self.assertFalse(result["legacy_method"]["available"])
        self.assertEqual(result["resolved_candidate_variants"], 0)

    def test_descriptive_panels_are_reported_without_association_or_causal_claims(self):
        lead = {**self.lead, "selected_markers": "3", "status": "analysed"}
        table(self.run / "lead_loci.tsv", list(lead), [lead])
        result = audit_finemapping(self.run, self.gwas, self.legacy, self.out)
        self.assertEqual(result["cohort_selected_marker_panels"], 1)
        self.assertEqual(result["cohort_haplotype_association_tests"], 0)
        self.assertEqual(result["resolved_candidate_variants"], 0)
        self.assertIn("1 selected marker panels (descriptive only)", (self.out / "REPORT.txt").read_text())

    def test_tested_table_cannot_be_hidden_by_zero_summary_count(self):
        table(self.run / "haplotype_comparisons.tsv", ["test_status"], [{"test_status": "tested_unadjusted"}])
        with self.assertRaisesRegex(ValueError, "tested haplotypes"):
            audit_finemapping(self.run, self.gwas, self.legacy, self.out)
        self.assertFalse(self.out.exists())

    def test_synthetic_incomplete_or_tested_cohorts_are_rejected(self):
        for update in ({"synthetic": True}, {"haplotype_carrier_tests_performed": 1},
                       {"haplotype_overall_tests_performed": 1}, {"inputs_stable_during_read": False}):
            (self.run / "run_summary.json").write_text(json.dumps({**self.summary, **update}))
            with self.assertRaises(ValueError):
                audit_finemapping(self.run, self.gwas, self.legacy, self.out)
        (self.run / "COMPLETED.txt").unlink()
        with self.assertRaisesRegex(ValueError, "incomplete"):
            audit_finemapping(self.run, self.gwas, self.legacy, self.out)
        self.assertFalse(self.out.exists())


if __name__ == "__main__":
    unittest.main()
