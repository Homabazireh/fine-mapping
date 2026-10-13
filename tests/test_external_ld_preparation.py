"""Synthetic checks for external LD: no cohort findings are produced."""
from pathlib import Path
import csv
import gzip
import json
import shutil
import subprocess
import tempfile
import unittest

import numpy as np

from prepare_external_finemapping import (
    MATCH_FIELDS, common_af, complete_dosages, match_variants, normalize_gwas,
    prepare, read_reference, read_tsv, sha, signed_ld, variant_key, write_tsv,
)


def gwas_row(pos, ref="A", alt="G", af=.5):
    return dict(locus_id="L0001", variant_id=f"chr1:{pos}:{ref}:{alt}", chrom="chr1", pos=str(pos),
                ref=ref, alt=alt, variant_type="SNP", beta="0.2", se="0.05", pval="0.8", af=str(af))


class GenotypeAndLDTests(unittest.TestCase):
    def test_observed_diploid_dosages_only(self):
        np.testing.assert_array_equal(complete_dosages(["0/0", "0|1", "1/0", "1|1"], 4), [0, 1, 1, 2])
        for bad in ["./.", "0/.", "0|2", "1", "0/1/1", "0:10"]:
            self.assertIsNone(complete_dosages(["0/0", bad], 2))
        with self.assertRaises(ValueError):
            complete_dosages(["0/0"], 2)

    def test_signed_correlation_not_squared_ld(self):
        matrix = signed_ld([[0, 1, 2, 0], [2, 1, 0, 2], [0, 0, 1, 2]])
        self.assertAlmostEqual(matrix[0, 1], -1)
        np.testing.assert_allclose(matrix, matrix.T)
        self.assertGreaterEqual(np.linalg.eigvalsh(matrix).min(), -1e-12)
        np.testing.assert_array_equal(np.diag(matrix), np.ones(3))
        for bad in [[[0, 0, 0], [0, 1, 2]], [[0, 1, np.nan], [0, 1, 2]]]:
            with self.assertRaises(ValueError):
                signed_ld(bad)

    def test_matching_preserves_effect_sign_and_order(self):
        rows = [gwas_row(20), gwas_row(10)]
        keys = {i: variant_key(r["chrom"], r["pos"], r["ref"], r["alt"]) for i, r in enumerate(rows)}
        ref = {keys[0]: (np.array([0, 1, 2]), []), keys[1]: (np.array([2, 1, 0]), [])}
        matched, x, qc = match_variants(rows, keys, ref, 1, 30, 100)
        self.assertEqual([r["pos"] for r in matched], [10, 20])
        self.assertEqual([r["beta"] for r in matched], ["0.2", "0.2"])
        np.testing.assert_array_equal(x[0], [2, 1, 0])
        self.assertTrue(all(r["status"] == "matched" for r in qc))
        self.assertEqual(len(matched), 2)  # p=.8 does not remove variants.

    def test_flip_missing_duplicate_and_frequency_blockers(self):
        rows = [gwas_row(10, "G", "A"), gwas_row(20), gwas_row(30), gwas_row(40, af=.005)]
        keys = {i: variant_key(r["chrom"], r["pos"], r["ref"], r["alt"]) for i, r in enumerate(rows)}
        ref = {variant_key("chr1", 10, "A", "G"): (np.array([0, 1, 2]), []),
               keys[1]: (None, ["reference_incomplete_or_nonbiallelic_GT"]),
               keys[2]: (None, ["reference_duplicate_normalized_key"]),
               keys[3]: (np.array([2, 1, 0]), [])}
        matched, dosages, qc = match_variants(rows, keys, ref, 1, 50, 100)
        self.assertFalse(matched)
        self.assertFalse(dosages)
        self.assertIn("no_exact_normalized_reference_match", qc[0]["reasons"])
        self.assertIn("reference_incomplete", qc[1]["reasons"])
        self.assertIn("reference_duplicate", qc[2]["reasons"])
        self.assertIn("GWAS_MAF", qc[3]["reasons"])
        self.assertIn("AF_difference", qc[3]["reasons"])

    def test_normalized_duplicate_excludes_every_source_record(self):
        rows = [gwas_row(10, "AA", "A"), gwas_row(20, "AA", "A")]
        key = variant_key("chr1", 1, "AA", "A")
        matched, _, qc = match_variants(rows, {0:key, 1:key}, {key:(np.array([0,1,2]), [])}, 1, 30, 100)
        self.assertFalse(matched)
        self.assertTrue(all("GWAS_duplicate_normalized_key" in r["reasons"] for r in qc))
        matched, _, qc = match_variants(rows[:1], {0:key}, {key:(np.array([0,1,2]), [])}, 2, 30, 100)
        self.assertFalse(matched)
        self.assertIn("outside_window", qc[0]["reasons"])

    def test_invalid_statistics_and_reference_monomorphism(self):
        row = gwas_row(10)
        row["af"] = "NA"
        key = variant_key("chr1", 10, "A", "G")
        matched, _, qc = match_variants([row], {0:key}, {key:(np.array([0,0,0]),[])}, 1, 20, 100)
        self.assertFalse(matched)
        self.assertIn("invalid_GWAS_summary_statistics", qc[0]["reasons"])
        self.assertIn("reference_MAF", qc[0]["reasons"])
        self.assertFalse(common_af("NA"))
        self.assertFalse(common_af(float("nan")))
        self.assertTrue(common_af(.01))

    def test_all_heterozygotes_excluded_despite_half_allele_frequency(self):
        row = gwas_row(10)
        key = variant_key("chr1", 10, "A", "G")
        matched, _, qc = match_variants([row], {0:key}, {key:(np.ones(503),[])}, 1, 20, 100)
        self.assertFalse(matched)
        self.assertEqual(qc[0]["reasons"], "reference_no_dosage_variation")


@unittest.skipUnless(shutil.which("bcftools"), "bcftools is required for normalization integration checks")
class BcftoolsIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.fasta = self.root / "reference.fa"
        self.fasta.write_text(">chr1\n" + "A"*100 + "\n")
        Path(str(self.fasta)+".fai").write_text("chr1\t100\t6\t100\t101\n")

    def tearDown(self):
        self.temp.cleanup()

    def test_indel_left_normalization_and_reference_mismatch_exclusion(self):
        rows = [gwas_row(20, "AA", "A"), gwas_row(30, "G", "A"), gwas_row(40, "A", "C")]
        normalized, command = normalize_gwas(rows, self.root, "chr1", self.fasta, 100)
        self.assertEqual(normalized[0], ("chr1", 1, "AA", "A"))
        self.assertNotIn(1, normalized)  # no automatic REF/ALT swap or beta flip.
        self.assertEqual(normalized[2], ("chr1", 40, "A", "C"))
        self.assertIn("x", command)

    def reference_fixture(self):
        reference = self.root/"reference"
        (reference/"chr1").mkdir(parents=True)
        samples = [f"S{i:04d}" for i in range(503)]
        (reference/"source_sample_panel.tsv").write_text("sample\tpop\tsuper_pop\tgender\n"+"".join(f"{s}\tGBR\tEUR\tmale\n" for s in samples))
        (reference/"EUR_503.samples.txt").write_text("\n".join(samples)+"\n")
        (reference/"chr1"/"regions.tsv").write_text("chr1\t1\t100\n")
        vcf = reference/"chr1"/"input.vcf"
        gt1 = ["0/0", "0/1", "1/1"]*167+["0/0", "0/1"]
        gt2 = ["1/1", "0/1", "0/0"]*167+["1/1", "0/1"]
        vcf.write_text("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"+
                      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+"\t".join(samples)+"\n"+
                      "chr1\t10\t.\tA\tG\t.\t.\t.\tGT\t"+"\t".join(gt1)+"\n"+
                      "chr1\t20\t.\tA\tC\t.\t.\t.\tGT\t"+"\t".join(gt2)+"\n"+
                      "chr1\t30\t.\tA\tT\t.\t.\t.\tGT\t"+"\t".join(gt1[:-1]+["./."])+"\n")
        normalized = reference/"chr1"/"normalized.vcf.gz"
        subprocess.run(["bcftools", "view", "-Oz", "-o", str(normalized), str(vcf)], check=True, capture_output=True)
        subprocess.run(["bcftools", "index", str(normalized)], check=True, capture_output=True)
        manifest = dict(genome_build="GRCh38", population="EUR_original_2504_unrelated_subset", sample_count=503,
            upstream_filtered_release_ACCEPTS_PASS_or_dot=True, fasta_sha256=sha(self.fasta),
            sample_panel_sha256=sha(reference/"source_sample_panel.tsv"), chromosomes=[dict(chrom="chr1",
                normalized_file="chr1/normalized.vcf.gz", normalized_sha256=sha(normalized),
                regions_sha256=sha(reference/"chr1"/"regions.tsv"))])
        (reference/"manifest.json").write_text(json.dumps(manifest))
        (reference/"COMPLETED.txt").write_text("Synthetic fixture only\n")
        return reference

    def test_full_preparation_preserves_alignment_and_excludes_missing_GT(self):
        reference = self.reference_fixture()
        source = self.root/"source"
        (source/"regions").mkdir(parents=True)
        (source/"lead_loci.tsv").write_text("locus_id\tchrom\tlead_pos\tstart\tend\nL0001\tchr1\t10\t1\t100\n")
        rows = [gwas_row(20, alt="C"), gwas_row(10), gwas_row(30, alt="T")]
        path = source/"regions"/"L0001.tsv"
        write_tsv(path, rows, list(rows[0]))
        manifest = dict(source_provenance=dict(source_build="GRCh38",effect_allele="ALT",
            case_count_from_portal=20,control_count_from_portal=80,total_participants_from_portal=100),
            source_files={"leads":{"sha256":sha(source/"lead_loci.tsv")}}, region_files={"L0001":{"sha256":sha(path)}})
        (source/"manifest.json").write_text(json.dumps(manifest))
        (source/"COMPLETED.txt").write_text("Synthetic fixture only\n")
        out = self.root/"prepared"
        result = prepare(source, reference, self.fasta, out, workers=1)
        self.assertEqual(result["prepared_regions"], 1)
        locus = result["loci"][0]
        variants = read_tsv(out/"L0001"/"matched_variants.tsv")
        self.assertEqual([r["variant_id"] for r in variants], locus["ld_variant_ids"])
        self.assertEqual([r["pos"] for r in variants], ["10", "20"])
        with gzip.open(out/"L0001"/"ld.tsv.gz", "rt") as stream:
            matrix = np.loadtxt(stream)
        self.assertEqual(matrix.shape, (2,2))
        self.assertAlmostEqual(matrix[0,1], -1)
        qc = read_tsv(out/"L0001"/"variant_qc.tsv")
        self.assertIn("reference_incomplete", qc[2]["reasons"])
        self.assertEqual(sha(out/"L0001"/"ld.tsv.gz"), locus["ld_sha256"])
        with self.assertRaises(ValueError):
            prepare(source, reference, self.fasta, out, workers=1)

    def test_reference_query_rejects_all_duplicate_keys_and_missing_GT(self):
        reference = self.reference_fixture()
        path = reference/"chr1"/"normalized.vcf.gz"
        original = reference/"chr1"/"input.vcf"
        lines = original.read_text().splitlines()
        duplicate = next(line for line in lines if line.startswith("chr1\t20\t"))
        lines.insert(lines.index(duplicate), duplicate)
        original.write_text("\n".join(lines)+"\n")
        subprocess.run(["bcftools", "view", "-Oz", "-o", str(path), str(original)], check=True, capture_output=True)
        subprocess.run(["bcftools", "index", "-f", str(path)], check=True, capture_output=True)
        records, _, counts = read_reference(path, "chr1", 1, 100, 503)
        self.assertEqual(len(records), 3)
        self.assertEqual(sum(counts.values()), 4)
        self.assertIsNone(records[("chr1",30,"A","T")][0])
        self.assertIsNone(records[("chr1",20,"A","C")][0])
        self.assertEqual(records[("chr1",20,"A","C")][1], ["reference_duplicate_normalized_key"])


if __name__ == "__main__":
    unittest.main()
