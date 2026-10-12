"""Verify exact-ID group denominators in aggregate lead-site genotype QC."""
from __future__ import annotations

import contextlib
import csv
import io
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'code'))

import check_lead_genotypes as checker


class LeadGenotypeTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='test_lead_genotypes_')
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)
        self.metadata = self.folder / 'categories.tsv'
        # Deliberately use a different order from the VCF. The similar-looking
        # unmatchedLR name must not classify the unmatched VCF sample.
        self.metadata.write_text(
            'name (Tubingen ID)\tCategory\n'
            'diabetes_high\tT2D\ncontrol_missing\tNORM\n'
            'diabetes_low\tT2D\nunmatchedLR\tPrediabetes\n'
        )
        self.record = (
            'chr1\t100\t.\tA\tG\t60\tPASS\t.\tGT:DP:GQ\t'
            '0/0:30:99\t0/0:30:5\t./.:.:.\t0/1:30:99\n'
        )
        vcf = self.folder / 'cohort.vcf'
        vcf.write_text(
            '##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
            'unmatched\tdiabetes_low\tcontrol_missing\tdiabetes_high\n' + self.record
        )
        gwas = self.folder / 'gwas.tsv'
        gwas.write_text('Chromsome\tPosition\tEffectAllele\tNonEffectAllele\tPval\n1\t100\tG\tA\t1e-12\n')
        self.config = self.folder / 'settings.json'
        self.config.write_text(json.dumps({'paths': {'vcf': str(vcf), 'gwas': str(gwas),
                                                   'metadata': str(self.metadata)}}))

    def run_check(self, out, metadata=None):
        query = SimpleNamespace(returncode=0, stdout=self.record, stderr='')
        with mock.patch.object(checker.shutil, 'which', return_value='/fake/bcftools'), \
                mock.patch.object(checker.subprocess, 'run', return_value=query), \
                contextlib.redirect_stdout(io.StringIO()):
            return checker.check_lead_genotypes(self.config, out, metadata_path=metadata)

    def test_group_counts_use_exact_sample_indices_and_keep_zero_distinct_from_na(self):
        output = self.folder / 'qc'
        summary = self.run_check(output, self.metadata)
        with (output / 'group_lead_genotype_qc.tsv').open() as handle:
            rows = {row['group']: row for row in csv.DictReader(handle, delimiter='\t')}
        self.assertEqual(set(rows), {'T2D', 'NORM', 'Prediabetes'})
        self.assertEqual(rows['T2D']['n_group'], '2')
        self.assertEqual(rows['T2D']['raw_called_samples'], '2')
        self.assertEqual(rows['T2D']['high_quality_called_samples'], '1')
        self.assertEqual(rows['T2D']['high_quality_call_rate'], '0.5')
        self.assertEqual(rows['T2D']['raw_reference_genotypes'], '1')
        self.assertEqual(rows['T2D']['high_quality_reference_genotypes'], '0')
        self.assertEqual(rows['NORM']['n_group'], '1')
        self.assertEqual(rows['NORM']['high_quality_call_rate'], '0')
        self.assertEqual(rows['Prediabetes']['n_group'], '0')
        self.assertEqual(rows['Prediabetes']['high_quality_called_samples'], '0')
        self.assertEqual(rows['Prediabetes']['high_quality_call_rate'], 'NA')
        self.assertEqual(summary['n_samples'], 4)
        self.assertEqual(summary['max_high_quality_called'], 2)
        self.assertEqual(summary['high_quality_reference_genotypes'], 1)
        group_qc = summary['group_qc']
        self.assertEqual(group_qc['matched_category_counts'], {'T2D': 2, 'NORM': 1, 'Prediabetes': 0})
        self.assertEqual(group_qc['source_category_counts'], {'T2D': 2, 'NORM': 1, 'Prediabetes': 1})
        self.assertEqual(group_qc['matched_sample_count'], 3)
        self.assertEqual(group_qc['unmatched_vcf_count'], 1)
        self.assertEqual(group_qc['unmatched_metadata_count'], 1)
        self.assertEqual(group_qc['scope'], 'exact_metadata_id_matches_only')
        self.assertEqual(Path(group_qc['table']), output / 'group_lead_genotype_qc.tsv')
        self.assertIn('sha256', summary['inputs']['metadata'])

    def test_metadata_is_explicit_and_does_not_change_all_sample_counts(self):
        plain = self.run_check(self.folder / 'plain')
        grouped = self.run_check(self.folder / 'grouped', self.metadata)
        self.assertNotIn('group_qc', plain)
        self.assertNotIn('metadata', plain['inputs'])
        self.assertFalse((self.folder / 'plain' / 'group_lead_genotype_qc.tsv').exists())
        self.assertEqual((self.folder / 'plain' / 'lead_genotype_qc.tsv').read_bytes(),
                         (self.folder / 'grouped' / 'lead_genotype_qc.tsv').read_bytes())
        for field in ['n_samples', 'n_leads', 'max_high_quality_called', 'raw_reference_genotypes',
                      'high_quality_reference_genotypes', 'lead_status_counts']:
            self.assertEqual(plain[field], grouped[field])

    def test_duplicate_metadata_cannot_produce_group_qc(self):
        with self.metadata.open('a') as handle:
            handle.write('diabetes_high\tNORM\n')
        output = self.folder / 'qc'
        with self.assertRaisesRegex(ValueError, 'duplicate metadata ID'):
            self.run_check(output, self.metadata)
        self.assertFalse(output.exists())

    def test_existing_preflight_outputs_are_preserved(self):
        output = self.folder / 'qc'
        output.mkdir()
        sentinel = output / 'summary.json'
        sentinel.write_text('{"previous": true}\n')
        with self.assertRaisesRegex(ValueError, 'Output directory is not empty'):
            self.run_check(output, self.metadata)
        self.assertEqual(sentinel.read_text(), '{"previous": true}\n')
        self.assertEqual(list(output.iterdir()), [sentinel])


if __name__ == '__main__':
    unittest.main()
