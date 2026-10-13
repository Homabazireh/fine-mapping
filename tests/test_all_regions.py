"""Every regional window remains visible without changing scientific denominators."""
from __future__ import annotations

import csv
from copy import deepcopy
import gzip
from html.parser import HTMLParser
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "code"))
from plot_all_regions import export_all_regions, load_run, prepare_plot_data, read_table, validate_variant


class LinkParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "a":
            self.links.append(attrs["href"])


class AllRegionsTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory(prefix="all_regions_")
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.run = self.root / "cohort"
        self.run.mkdir()
        self.out = self.root / "regional_figures"
        self.summary = {"synthetic": False, "inputs_stable_during_read": True,
            "group_counts": {"T2D": 2, "Prediabetes": 2, "Control": 4},
            "selected_loci": 3, "distance_clusters_available": 3, "genome_build_user_declared": "GRCh38",
            "parameters": {"top_loci": 0, "min_dp": 10, "min_gq": 20, "min_call_rate": .9}}
        self.leads = [{"locus_id": f"L{i:04d}", "chrom": "chr1", "lead_pos": str(100*i),
                       "start": str(100*i-50), "end": str(100*i+50), "status": "lead_below_min_call_rate"}
                      for i in range(1, 4)]
        (self.run / "COMPLETED.txt").write_text("Completed\n")
        self.rows = self.variant(self.leads[0], 90) + self.variant(self.leads[0], 100, exact=True)
        # The second region has only a non-lead SNP; the third has no biallelic SNP.
        self.rows += self.variant(self.leads[1], 210)
        self.rows += self.variant(self.leads[2], 300, alt="GA", kind="INDEL")
        self.save_run()

    def variant(self, lead, pos, *, exact=False, alt="G", kind="SNP"):
        result = []
        counts = {"T2D": (2, [3, 1]), "Prediabetes": (1, [1, 1]), "Control": (0, [0, 0])}
        for group, (called, copies) in counts.items():
            for index, allele in enumerate(["A", alt]):
                result.append({"locus_id": lead["locus_id"], "variant_id": f"chr1:{pos}:A:{alt}",
                    "chrom": "chr1", "pos": str(pos), "ref": "A", "alt": alt, "type": kind,
                    "distance_to_lead_bp": str(pos-int(lead["lead_pos"])), "is_exact_lead_match": str(int(exact)),
                    "group": group, "n_group": str(self.summary["group_counts"][group]), "n_called": str(called),
                    "called_chromosomes": str(2*called), "allele_index": str(index), "allele": allele,
                    "allele_copies": str(copies[index]), "call_rate": "0.987", "allele_frequency": "0.987"})
        return result

    def save_run(self):
        (self.run / "run_summary.json").write_text(json.dumps(self.summary))
        with (self.run / "lead_loci.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, list(self.leads[0]), delimiter="\t")
            writer.writeheader(); writer.writerows(self.leads)
        with gzip.open(self.run / "regional_variant_frequencies.tsv.gz", "wt", newline="") as handle:
            writer = csv.DictWriter(handle, list(self.rows[0]), delimiter="\t")
            writer.writeheader(); writer.writerows(self.rows)

    def test_partial_cluster_selection_and_unstable_inputs_are_rejected_before_write(self):
        for update in [{"distance_clusters_available": 5}, {"inputs_stable_during_read": False}]:
            with self.subTest(update=update):
                (self.run / "run_summary.json").write_text(json.dumps({**self.summary, **update}))
                with self.assertRaises(ValueError):
                    export_all_regions(self.run, self.out)
                self.assertFalse(self.out.exists())

    def test_synthetic_run_requires_explicit_fixture_opt_in(self):
        self.summary["synthetic"] = True
        self.save_run()
        with self.assertRaisesRegex(ValueError, "synthetic"):
            load_run(self.run)
        self.assertTrue(load_run(self.run, allow_synthetic=True)[0]["synthetic"])

    def test_nonlead_sites_and_zero_snp_regions_are_preserved(self):
        self.out.mkdir()
        metadata = prepare_plot_data(self.run, self.out, self.summary, self.leads)
        self.assertEqual([r["biallelic_snps"] for r in metadata], [2, 1, 0])
        self.assertEqual(metadata[2]["omitted_record_types"], {"INDEL": 1})
        second = read_table(self.out / metadata[1]["folder"] / "plot_data.tsv.gz")
        self.assertEqual({r["pos"] for r in second}, {"210"})
        self.assertEqual({r["is_exact_lead_match"] for r in second}, {"0"})
        self.assertEqual(read_table(self.out / metadata[2]["folder"] / "plot_data.tsv.gz"), [])

    def test_count_denominators_override_cached_fractions_and_zero_calls_are_na(self):
        rows = validate_variant(self.rows[:6], self.leads[0], self.summary["group_counts"])
        cells = {r["group"]: r for r in rows}
        self.assertEqual(cells["T2D"]["completeness"], 1)
        self.assertEqual(cells["T2D"]["observed_alt_fraction"], .25)
        self.assertEqual(cells["Control"]["completeness"], 0)
        self.assertEqual(cells["Control"]["observed_alt_fraction"], "NA")

    def test_bad_allele_sums_group_counts_or_distances_are_rejected(self):
        for field, value in [("allele_copies", "2"), ("called_chromosomes", "3"),
                             ("n_group", "200"), ("distance_to_lead_bp", "999")]:
            with self.subTest(field=field):
                rows = deepcopy(self.rows[:6]); rows[0][field] = value
                with self.assertRaises(ValueError):
                    validate_variant(rows, self.leads[0], self.summary["group_counts"])
        with self.assertRaisesRegex(ValueError, "group"):
            validate_variant(self.rows[:4], self.leads[0], self.summary["group_counts"])
        misplaced_lead = [{**r, "is_exact_lead_match": "1"} for r in self.rows[:6]]
        with self.assertRaisesRegex(ValueError, "lead flag"):
            validate_variant(misplaced_lead, self.leads[0], self.summary["group_counts"])

    def test_duplicate_blocks_cannot_silently_double_marker_counts(self):
        self.rows += self.rows[:6]
        self.save_run()
        self.out.mkdir()
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            prepare_plot_data(self.run, self.out, self.summary, self.leads)

    def test_indel_opt_in_preserves_missingness_and_snp_counts(self):
        self.out.mkdir()
        metadata = prepare_plot_data(self.run, self.out, self.summary, self.leads, include_indels=True)
        self.assertEqual([r["biallelic_snps"] for r in metadata], [2, 1, 0])
        self.assertEqual([r["biallelic_indels"] for r in metadata], [0, 0, 1])
        self.assertEqual([r["plotted_variant_count"] for r in metadata], [2, 1, 1])
        indels = read_table(self.out / metadata[2]["folder"] / "plot_data.tsv.gz")
        self.assertEqual({r["variant_type"] for r in indels}, {"indel"})
        control = next(r for r in indels if r["group"] == "Control")
        self.assertEqual(control["completeness"], "0.0")
        self.assertEqual(control["observed_alt_fraction"], "NA")
        self.assertEqual(metadata[2]["other_records_omitted"], 0)

    def test_indel_sequence_length_and_allele_conservation_are_checked(self):
        for ref, alt in [("A", "AT"), ("AT", "A"), ("A", "A" * 50)]:
            with self.subTest(ref=ref, alt=alt):
                rows = self.variant(self.leads[0], 100, exact=True, alt=alt, kind="indel")
                for row in rows:
                    row["ref"] = ref
                    row["allele"] = ref if row["allele_index"] == "0" else alt
                exported = validate_variant(rows, self.leads[0], self.summary["group_counts"], include_indels=True)
                self.assertEqual(len(exported), 3)
                self.assertEqual({r["variant_type"] for r in exported}, {"indel"})
                self.assertEqual(exported[0]["observed_alt_fraction"], .25)
                rows[0]["allele_copies"] = "99"
                with self.assertRaises(ValueError):
                    validate_variant(rows, self.leads[0], self.summary["group_counts"], include_indels=True)
        for alt in ["A" * 51, "<DEL>", "AN", "G"]:
            rows = self.variant(self.leads[0], 100, alt=alt, kind="indel")
            self.assertEqual(validate_variant(rows, self.leads[0], self.summary["group_counts"], include_indels=True), [])
        rows = self.variant(self.leads[0], 100, alt="AAA", kind="indel")
        self.assertEqual(validate_variant(rows, self.leads[0], self.summary["group_counts"],
                                         include_indels=True, max_indel_length=1), [])

    def test_mixed_multiallelic_record_remains_omitted_with_indel_opt_in(self):
        rows = self.variant(self.leads[0], 100)
        for row in rows:
            row["alt"] = "G,AT"
            row["type"] = "SNP/indel"
        for row in list(rows):
            if row["allele_index"] == "0":
                rows.append({**row, "allele_index": "2", "allele": "AT", "allele_copies": "0"})
        self.assertEqual(validate_variant(rows, self.leads[0], self.summary["group_counts"], include_indels=True), [])

    def test_indel_export_labels_classes_and_renders_indel_only_region(self):
        originals = {p.name: p.read_bytes() for p in self.run.iterdir()}
        manifest = export_all_regions(self.run, self.out, workers=1, dpi=72, include_indels=True)
        self.assertEqual(manifest["plotted_variant_count"], 4)
        self.assertEqual(manifest["biallelic_snps"], 3)
        self.assertEqual(manifest["biallelic_indels"], 1)
        region = manifest["regions"][2]
        text = (self.out / region["folder"] / "INTERPRETATION.md").read_text()
        self.assertIn("0 SNPs and 1 small indels", text)
        self.assertIn("triangles are small indels", text)
        self.assertIn("SNPs / indels", (self.out / "index.html").read_text())
        interpretation = read_table(self.out / "interpretations.tsv")
        third = [r for r in interpretation if r["locus_id"] == "L0003"]
        self.assertTrue(all(r["biallelic_indels"] == "1" and r["plotted_variant_count"] == "1" for r in third))
        self.assertEqual({p.name: p.read_bytes() for p in self.run.iterdir()}, originals)

    def test_separate_figures_and_interpretations_for_every_region_are_portable(self):
        try:
            import matplotlib  # noqa: F401
        except ImportError:
            self.skipTest("matplotlib not installed")
        originals = {p.name: p.read_bytes() for p in self.run.iterdir()}
        manifest = export_all_regions(self.run, self.out, workers=1, dpi=72)
        self.assertEqual(manifest["n_regions"], 3)
        self.assertEqual(manifest["n_figures"], 3)
        interpretation = read_table(self.out / "interpretations.tsv")
        self.assertEqual(len(interpretation), 9)
        for region in manifest["regions"]:
            folder = self.out / region["folder"]
            self.assertTrue((folder / "region.png").is_file())
            self.assertTrue((folder / "region.pdf").is_file())
            text = (folder / "INTERPRETATION.md").read_text()
            self.assertIn("may be biased", text)
            self.assertIn("do not show phased haplotypes, test group differences", text)
        blank = [r for r in interpretation if r["locus_id"] == "L0003"]
        self.assertTrue(all(r["median_completeness"] == "NA" for r in blank))
        self.assertEqual({p.name: p.read_bytes() for p in self.run.iterdir()}, originals)
        self.assertEqual((self.out / "all_regional_variant_frequencies.tsv.gz").read_bytes(),
                         (self.run / "regional_variant_frequencies.tsv.gz").read_bytes())
        html = (self.out / "index.html").read_text()
        parser = LinkParser(); parser.feed(html)
        relocated = self.root / "downloaded regional comparisons"
        self.out.rename(relocated)
        for link in parser.links:
            self.assertFalse(Path(link).is_absolute())
            self.assertTrue((relocated / link).is_file())
        self.assertNotIn(str(self.root), (relocated / "manifest.json").read_text())
        with self.assertRaisesRegex(ValueError, "new or empty"):
            export_all_regions(self.run, relocated)


if __name__ == "__main__":
    unittest.main()
