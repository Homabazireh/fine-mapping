"""Reference-download provenance checks with synthetic files and no network."""
import io
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
import fetch_external_ld_reference as fetch
from extract_external_gwas import fingerprint


class ReferenceFetchTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.leads = self.root / "lead_loci.tsv"
        self.leads.write_text("locus_id\tchrom\tstart\tend\tlead_pos\nL1\tchr1\t1\t4\t2\n")
        self.fasta = self.root / "reference.fa"
        self.fasta.write_text(">chr1\nACGT\n")
        self.index = Path(str(self.fasta) + ".fai")
        self.index.write_text("chr1\t4\t6\t4\t5\n")
        self.samples = [f"SYNTHETIC{i:03}" for i in range(503)]
        self.panel = ("sample\tpop\tsuper_pop\tgender\n" +
                      "".join(f"{s}\tTEST\tEUR\tunknown\n" for s in self.samples)).encode()

    def fake_run(self, command, **kwargs):
        if "-o" in command:
            Path(command[command.index("-o") + 1]).write_bytes(b"synthetic reference artifact")
        if command[1] == "index":
            Path(command[2] + ".csi").write_bytes(b"synthetic index")

    def fake_query(self, command, **kwargs):
        if command[1] == "query":
            return "\n".join(self.samples) + "\n"
        return "bcftools TEST\nUsing htslib TEST\n"

    def test_success_records_before_and_after_input_identity(self):
        out = self.root / "reference"
        with patch.object(fetch.urllib.request, "urlopen", return_value=io.BytesIO(self.panel)), \
             patch.object(fetch.subprocess, "run", side_effect=self.fake_run), \
             patch.object(fetch.subprocess, "check_output", side_effect=self.fake_query), \
             patch("sys.stdout", new_callable=io.StringIO):
            result = fetch.acquire(self.leads, self.fasta, out, workers=1)
        self.assertTrue((out / "COMPLETED.txt").exists())
        self.assertTrue(result["inputs_unchanged_during_acquisition"])
        self.assertEqual(result["fasta_sha256"], fingerprint(self.fasta)["sha256"])
        self.assertIn("fasta_index", result["input_snapshots"])

    def test_content_mutation_with_restored_metadata_blocks_completion(self):
        out = self.root / "mutated"
        before = self.fasta.stat()
        changed = False
        def mutate(command, **kwargs):
            nonlocal changed
            self.fake_run(command, **kwargs)
            if not changed:
                self.fasta.write_text(">chr1\nTGCA\n")
                os.utime(self.fasta, ns=(before.st_atime_ns, before.st_mtime_ns))
                changed = True
        with patch.object(fetch.urllib.request, "urlopen", return_value=io.BytesIO(self.panel)), \
             patch.object(fetch.subprocess, "run", side_effect=mutate), \
             patch.object(fetch.subprocess, "check_output", side_effect=self.fake_query), \
             patch("sys.stdout", new_callable=io.StringIO):
            with self.assertRaisesRegex(ValueError, "Input changed"):
                fetch.acquire(self.leads, self.fasta, out, workers=1)
        self.assertFalse((out / "COMPLETED.txt").exists())
        self.assertFalse((out / "manifest.json").exists())

    def test_changed_index_is_part_of_provenance_guard(self):
        snapshot = {"index": fingerprint(self.index)}
        self.index.write_text("chr1\t5\t6\t5\t6\n")
        with self.assertRaisesRegex(ValueError, "Input changed"):
            fetch.verify_input_snapshots(snapshot)

    def test_invalid_windows_and_worker_counts_fail_before_network(self):
        with patch.object(fetch.urllib.request, "urlopen") as network:
            for workers in [0, 9, 1.5]:
                with self.assertRaisesRegex(ValueError, "workers"):
                    fetch.acquire(self.leads, self.fasta, self.root / "bad", workers=workers)
            self.leads.write_text("locus_id\tchrom\tstart\tend\tlead_pos\nL1\tchr1\t1\t4\t9\n")
            with self.assertRaisesRegex(ValueError, "contain"):
                fetch.acquire(self.leads, self.fasta, self.root / "bad")
            network.assert_not_called()


if __name__ == "__main__":
    unittest.main()
