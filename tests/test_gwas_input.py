"""Regression checks for the DIAGRAM input and explicit synthetic GWAS overrides."""
from __future__ import annotations

import contextlib
import csv
from decimal import Decimal
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parent.parent
CODE = ROOT / "code"
sys.path.insert(0, str(CODE))

import compare_haplotypes as pipeline


REAL_GWAS = ROOT / "data" / "raw" / "EUR_diagram_t2d_1e-8.tsv"
DIAGRAM_HEADER = [
    "Chromsome", "Position", "EffectAllele", "NonEffectAllele",
    "Beta", "SE", "EAF", "Pval", "Ncases", "Ncontrols", "Neff", "liftover_strand",
]


class GwasInputTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="test_gwas_input_")
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)

    def arguments(self, *extra):
        return pipeline.argument_parser().parse_args([
            "--vcf", str(self.folder / "cohort.vcf"),
            "--samples", str(self.folder / "samples.tsv"),
            "--out", str(self.folder / "results"),
            "--genome-build", "GRCh38", *extra,
        ])

    def diagram_file(self, positions_and_p, chromosome="1"):
        path = self.folder / "diagram.tsv"
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(DIAGRAM_HEADER)
            for position, p_value in positions_and_p:
                writer.writerow([
                    chromosome, position, "G", "A", "0.2", "0.03", "0.4",
                    p_value, "100", "200", "266.67", "+",
                ])
        return path

    def test_diagram_schema_and_inclusive_default_cutoff(self):
        path = self.diagram_file([(100, "1e-12"), (110, "1e-8"), (120, "5e-8")])
        rows, counts, columns = pipeline.read_gwas(
            self.arguments("--gwas", str(path)), {"chr1": 1000},
        )
        self.assertEqual(columns, {
            "chrom": "Chromsome", "pos": "Position", "p": "Pval",
            "ea": "EffectAllele", "oa": "NonEffectAllele",
        })
        self.assertEqual([(row.chrom, row.pos) for row in rows], [("chr1", 100), ("chr1", 110)])
        self.assertEqual(counts["input_rows"], 3)
        self.assertEqual(counts["eligible_significant_rows"], 2)
        self.assertEqual(counts["above_p_threshold"], 1)

    def test_chromosome_translation_ignores_alternative_contigs(self):
        path = self.diagram_file([(100, "1e-12")])
        contigs = {f"chr1_KI{index}_alt": None for index in range(3000)}
        contigs.update({"chr1": 1000, "chrUn_KI270442v1": None, "chrM": None})
        rows, counts, _ = pipeline.read_gwas(self.arguments("--gwas", str(path)), contigs)
        self.assertEqual([row.chrom for row in rows], ["chr1"])
        self.assertEqual(counts["rows_with_contig_name_translation"], 1)

    def test_ambiguous_chromosome_aliases_are_rejected(self):
        path = self.diagram_file([(100, "1e-12")], chromosome="chr1")
        with self.assertRaisesRegex(ValueError, "Cannot map GWAS contig 'chr1' uniquely"):
            pipeline.read_gwas(self.arguments("--gwas", str(path)), {"1": None, "01": None})

    def test_exact_chromosome_match_takes_precedence_over_aliases(self):
        path = self.diagram_file([(100, "1e-12")])
        rows, counts, _ = pipeline.read_gwas(
            self.arguments("--gwas", str(path)), {"1": 1000, "chr1": 1000, "01": 1000},
        )
        self.assertEqual([row.chrom for row in rows], ["1"])
        self.assertEqual(counts.get("rows_with_contig_name_translation", 0), 0)

    def test_gwas_autosome_absent_from_vcf_is_rejected(self):
        path = self.diagram_file([(100, "1e-12")])
        with self.assertRaisesRegex(ValueError, "Cannot map GWAS contig '1' uniquely"):
            pipeline.read_gwas(self.arguments("--gwas", str(path)), {"chr2": None, "chrM": None})

    def test_explicit_threshold_override(self):
        path = self.diagram_file([(100, "1e-8"), (110, "5e-8")])
        args = self.arguments("--gwas", str(path), "--p-threshold", "5e-8")
        rows, counts, _ = pipeline.read_gwas(args, {"1": None})
        self.assertEqual([row.pos for row in rows], [100, 110])
        self.assertEqual(counts["eligible_significant_rows"], 2)
        self.assertEqual(Path(args.gwas), path)

    def test_legacy_column_aliases_with_explicit_input(self):
        path = self.folder / "legacy.tsv"
        path.write_text("CHR\tBP\tEA\tNEA\tP\n1\t100\tG\tA\t2e-8\n")
        rows, _, columns = pipeline.read_gwas(
            self.arguments("--gwas", str(path), "--p-threshold", "5e-8"), {"1": None},
        )
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].p, Decimal("2e-8"))
        self.assertEqual(columns["chrom"], "CHR")
        self.assertEqual(columns["p"], "P")

    def test_extremely_small_p_values_keep_decimal_ranking(self):
        path = self.diagram_file([(100, "1e-400"), (110, "3.35e-1003"), (120, "4e-1003")])
        rows, _, _ = pipeline.read_gwas(self.arguments("--gwas", str(path)), {"1": None})
        leads = pipeline.choose_leads(rows, gap=150000)
        self.assertEqual(len(leads), 1)
        self.assertEqual(leads[0][0].pos, 110)
        self.assertEqual(leads[0][0].p, Decimal("3.35e-1003"))
        self.assertTrue(all(row.p > 0 for row in rows))

    def test_default_gwas_is_independent_of_working_directory(self):
        previous = Path.cwd()
        try:
            os.chdir(self.folder)
            args = self.arguments()
        finally:
            os.chdir(previous)
        self.assertEqual(Path(args.gwas), REAL_GWAS)
        self.assertTrue(Path(args.gwas).is_absolute())
        self.assertEqual(args.p_threshold, "1e-8")

    def test_synthetic_end_to_end_uses_explicit_gwas(self):
        path = self.diagram_file([(100, "1e-12"), (110, "5e-8")])
        (self.folder / "samples.tsv").write_text(
            "sample\tgroup\nA1\tgroup_a\nA2\tgroup_a\nB1\tgroup_b\nB2\tgroup_b\n"
        )
        header = (
            "##fileformat=VCFv4.2\n##reference=GRCh38\n##contig=<ID=1,length=1000>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA1\tA2\tB1\tB2\n"
        )
        calls = "\t".join(["0|0:20:99:100", "0|1:20:99:100", "1|1:20:99:100", "0|1:20:99:100"])
        (self.folder / "cohort.vcf").write_text(
            header
            + f"1\t100\t.\tA\tG\t.\tPASS\t.\tGT:DP:GQ:PS\t{calls}\n"
            + f"1\t110\t.\tC\tT\t.\tPASS\t.\tGT:DP:GQ:PS\t{calls}\n"
        )
        args = self.arguments("--gwas", str(path))
        with contextlib.redirect_stdout(io.StringIO()):
            summary = pipeline.run(args)
        out = self.folder / "results"
        self.assertTrue((out / "COMPLETED.txt").is_file())
        self.assertEqual(summary["gwas_counts"]["input_rows"], 2)
        self.assertEqual(summary["gwas_counts"]["eligible_significant_rows"], 1)
        self.assertEqual(summary["locus_status_counts"], {"analysed": 1})
        self.assertEqual(summary["group_counts"], {"group_a": 2, "group_b": 2})
        self.assertEqual(summary["haplotype_carrier_tests_performed"], 0)
        saved = json.loads((out / "run_summary.json").read_text())
        self.assertEqual(saved["parameters"]["gwas"], str(path))
        with (out / "lead_loci.tsv").open() as handle:
            loci = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(loci), 1)
        self.assertEqual(loci[0]["lead_pos"], "100")
        self.assertEqual(loci[0]["selected_markers"], "2")

    @unittest.skipUnless(REAL_GWAS.is_file(), "Local DIAGRAM data is not distributed with the code")
    def test_local_diagram_file_counts_and_chromosome_translation(self):
        rows, counts, columns = pipeline.read_gwas(
            self.arguments(), {f"chr{chrom}": None for chrom in range(1, 23)},
        )
        self.assertGreater(counts["input_rows"], 0)
        self.assertEqual(counts["eligible_significant_rows"], len(rows))
        self.assertGreater(len(rows), 0)
        self.assertGreaterEqual(counts["input_rows"], len(rows) + counts.get("above_p_threshold", 0))
        self.assertGreater(counts.get("rows_with_contig_name_translation", 0), 0)
        self.assertLessEqual(counts["rows_with_contig_name_translation"], counts["input_rows"])
        self.assertEqual(columns["chrom"], "Chromsome")
        threshold = Decimal(self.arguments().p_threshold)
        self.assertTrue(all(0 <= row.p <= threshold for row in rows))
        self.assertTrue(all(row.chrom.startswith("chr") for row in rows))


if __name__ == "__main__":
    unittest.main()
