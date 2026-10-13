"""Synthetic external statistics only; no cohort files or network access."""
import csv
import gzip
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from extract_external_gwas import PROFILE, SOURCE_FIELDS, classify, extract, read_windows, valid_statistics


class ExternalExtractionTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.leads = self.root / "lead_loci.tsv"
        self.gwas = self.root / "gwas.tsv"
        self.out = self.root / "extracted"
        self.windows = [dict(locus_id="L1", chrom="chr1", lead_pos=100, start=90, end=110),
                        dict(locus_id="EMPTY", chrom="chr2", lead_pos=100, start=90, end=110)]
        self.write_windows()

    def write_windows(self):
        with self.leads.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, list(self.windows[0]), delimiter="\t")
            writer.writeheader(); writer.writerows(self.windows)

    def row(self, pos=100, ref="A", alt="G", **changes):
        row = dict(chrom="1", pos=str(pos), ref=ref, alt=alt, beta="0.123456789012345678",
                   sebeta="0.1", pval="0.95", af="0.2", case_af="0.3", control_af="0.19",
                   rsids="", nearest_genes="GENE", consequence="intron")
        return {**row, **changes}

    def write_gwas(self, rows, compressed=False):
        opener = gzip.open if compressed else open
        with opener(self.gwas, "wt", newline="") as stream:
            writer = csv.DictWriter(stream, sorted(SOURCE_FIELDS), delimiter="\t")
            writer.writeheader(); writer.writerows(rows)

    def extracted_rows(self):
        with gzip.open(self.out / "regional_variants.tsv.gz", "rt") as stream:
            return list(csv.DictReader(stream, delimiter="\t"))

    def test_inclusive_windows_snps_indels_precision_and_empty_regions(self):
        self.write_gwas([self.row(89), self.row(90), self.row(100, alt="AT", pval="1e-2000"),
                         self.row(110, ref="AT", alt="A"), self.row(111)], compressed=True)
        manifest = extract(self.gwas, self.leads, self.out)
        rows = self.extracted_rows()
        self.assertEqual([int(r["pos"]) for r in rows], [90, 100, 110])
        self.assertEqual([r["variant_type"] for r in rows], ["SNP", "small_indel", "small_indel"])
        self.assertEqual(rows[1]["pval"], "1e-2000")
        self.assertEqual(rows[0]["pval"], "0.95")  # No significance filter.
        self.assertEqual(rows[0]["beta"], "0.123456789012345678")
        self.assertEqual((self.out / "lead_loci.tsv").read_bytes(), self.leads.read_bytes())
        self.assertEqual((self.out / "regions/EMPTY.tsv").read_text().count("\n"), 1)
        self.assertTrue((self.out / "COMPLETED.txt").exists())
        self.assertFalse(manifest["source_provenance"]["byte_identity_to_remote_release_verified"])
        self.assertFalse(manifest["p_value_filter_applied"])
        self.assertEqual(manifest["new_statistical_tests_or_fits"], 0)
        self.assertIn("lead_loci.tsv", manifest["output_files"])

    def test_plain_text_and_distinct_alleles_at_same_position_are_preserved(self):
        self.write_gwas([self.row(100, alt="G"), self.row(100, alt="T")])
        extract(self.gwas, self.leads, self.out, source=PROFILE)
        self.assertEqual(len({r["variant_id"] for r in self.extracted_rows()}), 2)

    def test_sequence_and_length_gates_exclude_unsupported_alleles(self):
        self.write_gwas([self.row(90, alt="A"*50), self.row(91, alt="A"*51),
                         self.row(92, ref="AT", alt="GC"), self.row(93, alt="G,T"),
                         self.row(94, alt="<DEL>"), self.row(95, alt="N"), self.row(96, alt="A")])
        manifest = extract(self.gwas, self.leads, self.out)
        self.assertEqual([r["pos"] for r in self.extracted_rows()], ["90"])
        self.assertEqual(manifest["extraction_counts"]["excluded_MNV"], 1)
        self.assertEqual(classify("AC", "G"), "small_indel")  # Stored unequal sequence, not asserted pure deletion.

    def test_invalid_statistics_and_tiny_negative_p_values_are_excluded(self):
        self.write_gwas([self.row(90, pval="-1e-2000"), self.row(91, sebeta="0"),
                         self.row(92, beta="NaN"), self.row(93, af="1.01"), self.row(94)])
        manifest = extract(self.gwas, self.leads, self.out)
        self.assertEqual([r["pos"] for r in self.extracted_rows()], ["94"])
        self.assertEqual(manifest["extraction_counts"]["excluded_invalid_summary_statistics"], 4)
        self.assertTrue(valid_statistics(self.row(pval="1e-2000")))

    def test_duplicate_exact_variants_abort_even_if_first_statistics_invalid(self):
        self.write_gwas([self.row(sebeta="0"), self.row(chrom="chr1")])
        with self.assertRaisesRegex(ValueError, "Duplicate exact"):
            extract(self.gwas, self.leads, self.out)
        self.assertFalse((self.out / "COMPLETED.txt").exists())

    def test_invalid_duplicate_unsafe_or_overlapping_windows_are_rejected_before_write(self):
        self.write_gwas([self.row()])
        for changed in [dict(locus_id="L1"), dict(locus_id="../escape"), dict(start=0),
                        dict(lead_pos=200), dict(chrom="chrX"), dict(chrom="chr1", start=110)]:
            with self.subTest(changed=changed):
                original = self.windows[1].copy()
                self.windows[1].update(changed); self.write_windows()
                with self.assertRaises(ValueError):
                    extract(self.gwas, self.leads, self.out)
                self.assertFalse(self.out.exists())
                self.windows[1] = original

    def test_touching_inclusive_windows_overlap_but_adjacent_windows_do_not(self):
        self.windows[1] = dict(locus_id="L2", chrom="chr1", start=110, end=130, lead_pos=120)
        self.write_windows()
        with self.assertRaisesRegex(ValueError, "nonoverlapping"):
            read_windows(self.leads)
        self.windows[1]["start"] = 111
        self.write_windows()
        self.assertEqual(len(read_windows(self.leads)[0]), 2)

    def test_invalid_header_or_truncated_rows_never_get_completion_marker(self):
        self.gwas.write_text("chrom\tpos\n1\t100\n")
        with self.assertRaisesRegex(ValueError, "header"):
            extract(self.gwas, self.leads, self.out)
        self.assertFalse(self.out.exists())
        self.write_gwas([self.row()])
        with self.gwas.open("a") as stream:
            stream.write("incomplete\trow\n")
        with self.assertRaisesRegex(ValueError, "Malformed"):
            extract(self.gwas, self.leads, self.out)
        self.assertFalse((self.out / "COMPLETED.txt").exists())

    def test_existing_outputs_bad_limits_and_wrong_profile_are_rejected(self):
        self.write_gwas([self.row()])
        for limit in [0, 50, 1.5]:
            with self.assertRaisesRegex(ValueError, "max_indel_length"):
                extract(self.gwas, self.leads, self.out, max_indel_length=limit)
        with self.assertRaisesRegex(ValueError, "profile"):
            extract(self.gwas, self.leads, self.out, source="unknown")
        self.out.mkdir(); (self.out / "keep").write_text("original")
        with self.assertRaisesRegex(ValueError, "new output"):
            extract(self.gwas, self.leads, self.out)
        self.assertEqual((self.out / "keep").read_text(), "original")


if __name__ == "__main__":
    unittest.main()
