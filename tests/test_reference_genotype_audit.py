"""Reference-call audits must distinguish true reference, partial and missing GT."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "code"))
from audit_reference_genotypes import classify_gt_counts, merge_windows, sanitized_header, summarize_lead_lines


class ReferenceGenotypeAuditTests(unittest.TestCase):
    def test_raw_reference_survives_without_quality_filter(self):
        result = classify_gt_counts({"0/0": 3, "0|0": 2, "0|1": 4, "1/2": 1,
                                     "2/2": 7, "./.": 6, "./2": 8, "0": 1})
        self.assertEqual(result["ref_homozygotes"], 5)
        self.assertEqual(result["heterozygotes"], 5)
        self.assertEqual(result["alt_homozygotes"], 7)
        self.assertEqual(result["fully_missing_GT"], 6)
        self.assertEqual(result["partially_missing_GT"], 8)
        self.assertEqual(result["non_diploid_GT"], 1)
        self.assertEqual(result["all_GT_cells"], 32)
        self.assertEqual(result["raw_called_GT"] + result["missing_GT"], 32)

    def test_invalid_genotype_rejected(self):
        for genotype in ("", "0/x", "-1/0"):
            with self.assertRaises(ValueError):
                classify_gt_counts({genotype: 1})

    def test_overlapping_windows_counted_once(self):
        rows = [{"chrom": c, "start": s, "end": e} for c, s, e in
                [("chr1", 10, 20), ("chr1", 19, 25), ("chr1", 26, 30), ("chr2", 19, 25)]]
        self.assertEqual(merge_windows(rows), [("chr1", 10, 30), ("chr2", 19, 25)])

    def test_header_never_exports_sample_identifiers(self):
        header = ("##source=Caller\n##SAMPLE=<ID=PRIVATE_ID>\n"
                  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPRIVATE_ID\n")
        provenance, n = sanitized_header(header)
        self.assertEqual(provenance, ["##source=Caller"])
        self.assertEqual(n, 1)
        self.assertNotIn("PRIVATE_ID", str(provenance))

    def test_raw_lead_reference_is_counted_even_when_DP_fails(self):
        leads = [{"chrom": "chr1", "lead_pos": "12", "lead_effect_allele": "A", "lead_other_allele": "G"}]
        lines = ["chr1\t12\tA\tG\tPASS\t0/0:2:10:.\t./.:.:.:.\t0|1:20:30:1\n",
                 "chr1\t12\tA\tG,T\tPASS\t1/1:20:30:.\t1/1:20:30:.\t1/1:20:30:.\n",
                 "chr1\t12\tA\tG\tLowQual\t1/1:20:30:.\t1/1:20:30:.\t1/1:20:30:.\n"]
        result = summarize_lead_lines(lines, leads, 10, 20, 3)
        self.assertEqual(result["counts"]["matched_records"], 1)
        self.assertEqual(result["raw_GT_counts"]["ref_homozygotes"], 1)
        self.assertEqual(result["counts"]["called_DP_missing_or_lt_threshold"], 1)
        self.assertEqual(result["counts"]["called_DP_GQ_pass"], 1)
        self.assertEqual(result["missing_field_counts"]["missing_GT_DP_absent_GQ_absent"], 1)


if __name__ == "__main__":
    unittest.main()
