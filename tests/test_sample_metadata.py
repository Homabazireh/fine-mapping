"""Check that supplied groups can only be assigned to verified sample IDs."""
from __future__ import annotations

import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parent.parent
CODE = ROOT / "code"
sys.path.insert(0, str(CODE))

import make_sample_template as template


class SampleMetadataTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="test_sample_metadata_")
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)

    def table(self, name, columns, rows):
        path = self.folder / name
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(columns)
            writer.writerows(rows)
        return path

    def metadata(self, rows):
        return self.table("metadata.tsv", ["name (Tubingen ID)", "Category"], rows)

    def aliases(self, rows):
        return self.table("aliases.tsv", ["metadata_sample", "vcf_sample"], rows)

    def vcf(self, samples):
        path = self.folder / "cohort.vcf"
        path.write_text(
            "##fileformat=VCFv4.2\n##contig=<ID=1,length=1000>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            + "\t".join(samples) + "\n"
        )
        return path

    def cli(self, *arguments):
        return subprocess.run(
            [sys.executable, str(CODE / "make_sample_template.py"), *map(str, arguments)],
            cwd=self.folder, text=True, capture_output=True, check=False,
        )

    def test_exact_ids_only_and_values_preserved_with_whitespace_headers(self):
        path = self.table(
            "metadata.tsv", [" name (Tubingen ID) ", " Category ", "Short/Long "],
            [["exact_sample", "Prediabetes", "Long"],
             ["donorLR_01", "T2D", "Long"],
             [" leading_space", "NORM", "Long"],
             ["category_sample", "Mixed Case Group", "Long"]],
        )
        groups = template.metadata_groups(
            path, ["donor_02", "leading_space", "exact_sample", "category_sample", "absent"],
        )
        self.assertEqual(groups, {"exact_sample": "Prediabetes", "category_sample": "Mixed Case Group"})

    def test_explicit_alias_supports_different_prefix_suffix_and_identity(self):
        metadata = self.metadata([["donorLR_01", "T2D"], ["exact_sample", "NORM"]])
        aliases = self.aliases([["donorLR_01", "donor_02"], ["exact_sample", "exact_sample"]])
        self.assertEqual(
            template.metadata_groups(metadata, ["donor_02", "exact_sample"], aliases),
            {"donor_02": "T2D", "exact_sample": "NORM"},
        )

    def test_final_table_schema_preserves_ids_and_group_labels(self):
        path = self.table(
            'final.tsv', [' External Name ', ' Diabetes_Group ', ' name '],
            [['EXTERNAL_DEMO_001', 'Prediabetes', 'demoLRa001_01'],
             ['EXTERNAL_DEMO_002', 'Diabetes', 'with_space '],
             ['EXTERNAL_DEMO_003', 'Control', 'control']],
        )
        self.assertEqual(template.read_metadata_groups(path), {
            'demoLRa001_01': 'Prediabetes', 'with_space ': 'Diabetes', 'control': 'Control',
        })
        self.assertEqual(template.metadata_groups(path, ['control', 'EXTERNAL_DEMO_001', 'with_space']),
                         {'control': 'Control'})

    def test_both_metadata_schemas_must_agree_exactly(self):
        columns = ['name', 'Diabetes_Group', 'name (Tubingen ID)', 'Category']
        valid = self.table('valid.tsv', columns, [['sample', 'Diabetes', 'sample', 'Diabetes']])
        self.assertEqual(template.read_metadata_groups(valid), {'sample': 'Diabetes'})
        for row in [['sample', 'Diabetes', 'other', 'Diabetes'],
                    ['sample', 'Diabetes', 'sample', 'T2D'],
                    ['sample', 'Diabetes', 'sample ', 'Diabetes']]:
            with self.subTest(row=row):
                path = self.table('invalid.tsv', columns, [row])
                with self.assertRaisesRegex(ValueError, 'conflicting'):
                    template.read_metadata_groups(path)

    def test_incomplete_metadata_schema_is_rejected_even_without_rows(self):
        for columns in [['name', 'Category'], ['External Name', 'Diabetes_Group']]:
            with self.subTest(columns=columns):
                with self.assertRaisesRegex(ValueError, 'missing metadata schema'):
                    template.read_metadata_groups(self.table('invalid.tsv', columns, []))

    def test_alias_status_must_confirm_the_match_when_present(self):
        metadata = self.metadata([['source', 'Diabetes']])
        for status in ['pending_user_confirmation', '', 'unverified']:
            with self.subTest(status=status):
                aliases = self.table('aliases.tsv', ['metadata_sample', 'vcf_sample', 'status'],
                                     [['source', 'target', status]])
                with self.assertRaisesRegex(ValueError, 'status must be confirmed or verified'):
                    template.metadata_groups(metadata, ['target'], aliases)
        for status in ['confirmed', 'Verified', ' confirmed ']:
            with self.subTest(status=status):
                aliases = self.table('aliases.tsv', ['metadata_sample', 'vcf_sample', 'status'],
                                     [['source', 'target', status]])
                self.assertEqual(template.metadata_groups(metadata, ['target'], aliases),
                                 {'target': 'Diabetes'})

    def test_duplicate_or_empty_metadata_is_rejected(self):
        invalid_rows = {
            "conflicting duplicate": [["donor", "NORM"], ["donor", "T2D"]],
            "same group duplicate": [["donor", "NORM"], ["donor", "NORM"]],
            "empty ID": [["", "NORM"]],
            "blank ID": [["  ", "NORM"]],
            "empty category": [["donor", ""]],
            "blank category": [["donor", "  "]],
        }
        for label, rows in invalid_rows.items():
            with self.subTest(label=label):
                with self.assertRaises(ValueError):
                    template.metadata_groups(self.metadata(rows), ["donor"])

    def test_aliases_cannot_reference_unknown_ids(self):
        metadata = self.metadata([["source", "T2D"]])
        for row in [["unknown_source", "target"], ["source", "unknown_target"], ["", "target"]]:
            with self.subTest(alias=row):
                with self.assertRaises(ValueError):
                    template.metadata_groups(metadata, ["target"], self.aliases([row]))

    def test_aliases_must_be_one_to_one_even_when_groups_agree(self):
        metadata = self.metadata([["source_a", "T2D"], ["source_b", "T2D"]])
        invalid_maps = [
            [["source_a", "target_a"], ["source_a", "target_b"]],
            [["source_a", "target_a"], ["source_b", "target_a"]],
            [["source_a", "target_a"], ["source_a", "target_a"]],
        ]
        for rows in invalid_maps:
            with self.subTest(aliases=rows):
                with self.assertRaises(ValueError):
                    template.metadata_groups(metadata, ["target_a", "target_b"], self.aliases(rows))

    def test_aliases_cannot_override_or_redirect_exact_matches(self):
        metadata = self.metadata([["exact", "NORM"], ["alias_source", "T2D"]])
        for row in [["alias_source", "exact"], ["exact", "other_vcf_sample"]]:
            with self.subTest(alias=row):
                with self.assertRaises(ValueError):
                    template.metadata_groups(metadata, ["exact", "other_vcf_sample"], self.aliases([row]))

    def test_cli_writes_verified_groups_unresolved_ids_and_refuses_overwrite(self):
        vcf = self.vcf(["unmatched", "donor_02", "exact"])
        metadata = self.metadata([["donorLR_01", "T2D"], ["exact", "Prediabetes"]])
        aliases = self.aliases([["donorLR_01", "donor_02"]])
        output = self.folder / "samples.tsv"
        args = ["--vcf", vcf, "--metadata", metadata, "--id-map", aliases, "--out", output]
        result = self.cli(*args)
        self.assertEqual(result.returncode, 0, result.stderr)
        with output.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(rows, [
            {"sample": "unmatched", "group": "UNRESOLVED"},
            {"sample": "donor_02", "group": "T2D"},
            {"sample": "exact", "group": "Prediabetes"},
        ])
        original = output.read_bytes()
        second = self.cli(*args)
        self.assertNotEqual(second.returncode, 0)
        self.assertIn("Refusing to overwrite", second.stderr)
        self.assertEqual(output.read_bytes(), original)

    def test_cli_requires_metadata_for_aliases(self):
        output = self.folder / "samples.tsv"
        result = self.cli(
            "--vcf", self.vcf(["target"]), "--id-map", self.aliases([["source", "target"]]),
            "--out", output,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("--id-map requires --metadata", result.stderr)
        self.assertFalse(output.exists())

    def test_cli_without_metadata_preserves_unresolved_template_behavior(self):
        output = self.folder / "samples.tsv"
        result = self.cli("--vcf", self.vcf(["sample_b", "sample_a"]), "--out", output)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(output.read_text(), "sample\tgroup\nsample_b\tUNRESOLVED\nsample_a\tUNRESOLVED\n")


if __name__ == "__main__":
    unittest.main()
