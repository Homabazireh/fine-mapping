"""Verify exact-ID group denominators in aggregate lead-site genotype QC."""
from __future__ import annotations

import contextlib
import csv
import hashlib
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

    def run_check(self, out, metadata=None, id_map=None):
        query = SimpleNamespace(returncode=0, stdout=self.record, stderr='')
        with mock.patch.object(checker.shutil, 'which', return_value='/fake/bcftools'), \
                mock.patch.object(checker.subprocess, 'run', return_value=query), \
                contextlib.redirect_stdout(io.StringIO()):
            return checker.check_lead_genotypes(self.config, out, metadata_path=metadata, id_map_path=id_map)

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
        self.assertEqual(group_qc['exact_match_category_counts'], group_qc['matched_category_counts'])
        self.assertEqual(group_qc['verified_alias_category_counts'], {'T2D': 0, 'NORM': 0, 'Prediabetes': 0})
        self.assertEqual(group_qc['exact_match_sample_count'], 3)
        self.assertEqual(group_qc['verified_alias_sample_count'], 0)
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

    def test_verified_alias_uses_actual_vcf_column_and_records_provenance(self):
        aliases = self.folder / 'aliases.tsv'
        aliases.write_text('metadata_sample\tvcf_sample\tstatus\n'
                           'unmatchedLR\tunmatched\tconfirmed\n'
                           'diabetes_high\tdiabetes_high\tverified\n')
        exact = self.run_check(self.folder / 'exact', self.metadata)
        output = self.folder / 'aliases'
        summary = self.run_check(output, self.metadata, aliases)
        with (output / 'group_lead_genotype_qc.tsv').open() as handle:
            rows = {row['group']: row for row in csv.DictReader(handle, delimiter='\t')}
        # unmatched is the first VCF column, while unmatchedLR is last in metadata.
        self.assertEqual(rows['Prediabetes']['n_group'], '1')
        self.assertEqual(rows['Prediabetes']['high_quality_called_samples'], '1')
        self.assertEqual(rows['Prediabetes']['high_quality_call_rate'], '1')
        self.assertEqual(rows['Prediabetes']['high_quality_reference_genotypes'], '1')
        group_qc = summary['group_qc']
        self.assertEqual(group_qc['scope'], 'exact_and_verified_id_matches')
        self.assertEqual(group_qc['exact_match_category_counts'], {'T2D': 2, 'NORM': 1, 'Prediabetes': 0})
        self.assertEqual(group_qc['verified_alias_category_counts'], {'T2D': 0, 'NORM': 0, 'Prediabetes': 1})
        self.assertEqual(group_qc['matched_category_counts'], {'T2D': 2, 'NORM': 1, 'Prediabetes': 1})
        self.assertEqual(group_qc['exact_match_sample_count'], 3)
        self.assertEqual(group_qc['verified_alias_sample_count'], 1)
        self.assertEqual(group_qc['matched_sample_count'], 4)
        self.assertEqual(group_qc['unmatched_vcf_count'], 0)
        self.assertEqual(group_qc['unmatched_metadata_count'], 0)
        info = summary['inputs']['id_map']
        self.assertEqual(info['path'], str(aliases.resolve()))
        self.assertEqual(info['sha256'], hashlib.sha256(aliases.read_bytes()).hexdigest())
        self.assertEqual(info['size_bytes'], aliases.stat().st_size)
        self.assertEqual(info['mtime_ns'], aliases.stat().st_mtime_ns)
        self.assertEqual((self.folder / 'exact' / 'lead_genotype_qc.tsv').read_bytes(),
                         (output / 'lead_genotype_qc.tsv').read_bytes())
        for field in ['n_samples', 'n_leads', 'max_high_quality_called', 'raw_reference_genotypes',
                      'high_quality_reference_genotypes', 'lead_status_counts']:
            self.assertEqual(exact[field], summary[field])

    def test_final_metadata_schema_is_supported_in_group_qc(self):
        self.metadata.write_text('External Name\tDiabetes_Group\tname\n'
                                 '101\tDiabetes\tdiabetes_high\n'
                                 '102\tPrediabetes\tunmatchedLR\n'
                                 '103\tControl\tcontrol_missing\n'
                                 '104\tDiabetes\tdiabetes_low\n')
        summary = self.run_check(self.folder / 'final', self.metadata)
        self.assertEqual(summary['group_qc']['source_category_counts'],
                         {'Diabetes': 2, 'Prediabetes': 1, 'Control': 1})
        self.assertEqual(summary['group_qc']['matched_category_counts'],
                         {'Diabetes': 2, 'Prediabetes': 0, 'Control': 1})

    def test_alias_requires_explicit_metadata(self):
        with self.assertRaisesRegex(ValueError, '--id-map requires --metadata'):
            self.run_check(self.folder / 'invalid', id_map=self.folder / 'absent.tsv')
        with contextlib.redirect_stderr(io.StringIO()) as stderr:
            with self.assertRaises(SystemExit):
                checker.main(['--id-map', str(self.folder / 'absent.tsv')])
        self.assertIn('--id-map requires --metadata', stderr.getvalue())

    def test_invalid_or_unconfirmed_aliases_cannot_produce_outputs(self):
        invalid = [
            [('unmatchedLR', 'unmatched', 'pending_user_confirmation')],
            [('unmatchedLR', 'unmatched', 'confirmed'), ('unmatchedLR', 'unmatched', 'confirmed')],
            [('unmatchedLR', 'diabetes_high', 'confirmed')],
            [('diabetes_high', 'unmatched', 'confirmed')],
            [('absent', 'unmatched', 'confirmed')],
        ]
        for number, rows in enumerate(invalid):
            with self.subTest(rows=rows):
                aliases = self.folder / 'aliases.tsv'
                aliases.write_text('metadata_sample\tvcf_sample\tstatus\n' +
                                   ''.join('\t'.join(row) + '\n' for row in rows))
                output = self.folder / f'invalid_{number}'
                with self.assertRaises(ValueError):
                    self.run_check(output, self.metadata, aliases)
                self.assertFalse(output.exists())

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
