"""Check configuration paths and explicit CLI choices after project organization."""
from __future__ import annotations

import contextlib
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
import project_config


class ProjectConfigTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="test_project_config_")
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)

    def config_file(self, content):
        path = self.folder / "settings.json"
        path.write_text(json.dumps(content))
        return path

    @contextlib.contextmanager
    def unrelated_working_directory(self):
        previous = Path.cwd()
        try:
            os.chdir(self.folder)
            yield
        finally:
            os.chdir(previous)

    def test_organized_defaults_are_project_relative_from_other_directory(self):
        with self.unrelated_working_directory():
            config = project_config.load_config()
            args = pipeline.argument_parser().parse_args(["--genome-build", "GRCh38"])
        self.assertEqual(project_config.PROJECT_ROOT, ROOT)
        self.assertEqual(project_config.DEFAULT_CONFIG, ROOT / "config" / "settings.json")
        self.assertEqual(Path(config["paths"]["gwas"]), ROOT / "data" / "raw" / "EUR_diagram_t2d_1e-8.tsv")
        self.assertEqual(Path(config["paths"]["vcf"]), ROOT / "data" / "raw" / "full_LR-WGS.vcf.gz")
        self.assertEqual(Path(args.gwas), Path(config["paths"]["gwas"]))
        self.assertEqual(Path(args.vcf), Path(config["paths"]["vcf"]))
        self.assertTrue(Path(args.samples).is_absolute())
        self.assertTrue(Path(args.out).is_absolute())

    def test_custom_config_applies_defaults_and_cli_options_take_precedence(self):
        absolute_vcf = self.folder / "cohort.vcf"
        config = self.config_file({
            "paths": {"gwas": "data/raw/custom.tsv", "vcf": str(absolute_vcf)},
            "analysis": {"p_threshold": "2e-8", "top_loci": 3, "plot": True},
        })
        common = ["--genome-build", "GRCh38", "--config", str(config)]
        with self.unrelated_working_directory():
            configured = pipeline.argument_parser().parse_args(common)
            overridden = pipeline.argument_parser().parse_args(common + [
                "--gwas", str(self.folder / "override.tsv"), "--p-threshold", "4e-8",
                "--top-loci", "7", "--no-plot",
            ])
        self.assertEqual(Path(configured.gwas), ROOT / "data" / "raw" / "custom.tsv")
        self.assertEqual(Path(configured.vcf), absolute_vcf)
        self.assertEqual(configured.p_threshold, "2e-8")
        self.assertEqual(configured.top_loci, 3)
        self.assertTrue(configured.plot)
        self.assertEqual(Path(overridden.gwas), self.folder / "override.tsv")
        self.assertEqual(overridden.p_threshold, "4e-8")
        self.assertEqual(overridden.top_loci, 7)
        self.assertFalse(overridden.plot)

    def test_genome_build_still_requires_explicit_cli_declaration(self):
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as caught:
                pipeline.argument_parser().parse_args([])
        self.assertEqual(caught.exception.code, 2)

    def test_unknown_config_keys_cannot_silently_change_analysis(self):
        for content in [{"analysiz": {}}, {"analysis": {"p_thresold": "1e-5"}}]:
            with self.subTest(config=content):
                with self.assertRaises(ValueError):
                    project_config.load_config(self.config_file(content))


if __name__ == "__main__":
    unittest.main()
