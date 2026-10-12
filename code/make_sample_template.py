#!/usr/bin/env python3
"""Write VCF sample IDs with supplied metadata groups or UNRESOLVED labels."""
import argparse
import csv
from pathlib import Path

from compare_haplotypes import DEFAULT_VCF, vcf_header, write_table
from project_config import ConfigArgumentParser


def _tsv_rows(path, required_columns, alternative_schemas=()):
    """Read a TSV, tolerating surrounding header whitespace only."""
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if reader.fieldnames is None:
            raise ValueError(f'{path}: missing TSV header')
        columns = [column.strip() for column in reader.fieldnames]
        if len(columns) != len(set(columns)):
            raise ValueError(f'{path}: duplicate column names')
        missing = set(required_columns) - set(columns)
        if missing:
            raise ValueError(f'{path}: missing columns: {", ".join(sorted(missing))}')
        if alternative_schemas and not any(set(schema) <= set(columns) for schema in alternative_schemas):
            expected = ' or '.join(' + '.join(schema) for schema in alternative_schemas)
            raise ValueError(f'{path}: missing metadata schema: expected {expected}')
        reader.fieldnames = columns
        for row in reader:
            line = reader.line_num
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f'{path}: line {line} does not match the TSV header')
            yield line, row


def read_metadata_groups(metadata_path):
    """Read either participant-table schema without changing IDs or categories."""
    groups = {}
    schemas = [('name (Tubingen ID)', 'Category'), ('name', 'Diabetes_Group')]
    for line, row in _tsv_rows(metadata_path, [], alternative_schemas=schemas):
        sample_fields = [field for field in ('name (Tubingen ID)', 'name') if field in row]
        category_fields = [field for field in ('Category', 'Diabetes_Group') if field in row]
        for label, fields in [('sample ID', sample_fields), ('category', category_fields)]:
            if len({row[field] for field in fields}) != 1:
                raise ValueError(f'{metadata_path}: conflicting {label} columns on line {line}')
        sample = row[sample_fields[0]]
        category = row[category_fields[0]]
        if not sample.strip():
            raise ValueError(f'{metadata_path}: empty Tubingen ID on line {line}')
        if not category.strip():
            raise ValueError(f'{metadata_path}: empty Category on line {line}')
        if sample in groups:
            raise ValueError(f'{metadata_path}: duplicate metadata ID {sample!r}')
        groups[sample] = category
    return groups


def _resolve_metadata_groups(metadata, vcf_samples, id_map_path=None):
    """Match exact IDs, then apply explicitly supplied one-to-one aliases."""
    vcf_samples = list(vcf_samples)
    known_samples = set(vcf_samples)
    if len(known_samples) != len(vcf_samples):
        raise ValueError('Duplicate VCF sample IDs')
    groups = {sample: metadata[sample] for sample in vcf_samples if sample in metadata}
    if id_map_path is None:
        return groups

    sources, targets = set(), set()
    for line, row in _tsv_rows(id_map_path, ['metadata_sample', 'vcf_sample']):
        if 'status' in row and row['status'].strip().lower() not in {'confirmed', 'verified'}:
            raise ValueError(f'{id_map_path}: alias status must be confirmed or verified on line {line}')
        source, target = row['metadata_sample'], row['vcf_sample']
        if not source.strip() or not target.strip():
            raise ValueError(f'{id_map_path}: empty ID on line {line}')
        if source in sources:
            raise ValueError(f'{id_map_path}: duplicate metadata_sample {source!r}')
        if target in targets:
            raise ValueError(f'{id_map_path}: duplicate vcf_sample {target!r}')
        if source not in metadata:
            raise ValueError(f'{id_map_path}: metadata_sample {source!r} is absent from metadata')
        if target not in known_samples:
            raise ValueError(f'{id_map_path}: vcf_sample {target!r} is absent from the VCF')
        if source in known_samples and source != target:
            raise ValueError(f'{id_map_path}: cannot redirect exact VCF match {source!r}')
        if target in groups and target != source:
            raise ValueError(f'{id_map_path}: target {target!r} collides with an exact metadata match')
        sources.add(source)
        targets.add(target)
        groups[target] = metadata[source]
    return groups


def metadata_groups(metadata_path, vcf_samples, id_map_path=None):
    """Return VCF ID -> Category for exact or explicitly mapped metadata IDs.

    IDs and category values are preserved as supplied. Unmatched samples are
    absent from the result; callers should label them UNRESOLVED.
    """
    metadata = read_metadata_groups(metadata_path)
    return _resolve_metadata_groups(metadata, vcf_samples, id_map_path)


def main(argv=None):
    parser = ConfigArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  path_defaults={'vcf': 'vcf'})
    parser.add_argument('--vcf', default=str(DEFAULT_VCF), help='Multisample VCF or VCF.gz')
    parser.add_argument('--out', required=True, help='New sample metadata TSV')
    parser.add_argument('--metadata', help='TSV containing name + Diabetes_Group or name (Tubingen ID) + Category columns')
    parser.add_argument('--id-map', help='Verified aliases TSV with metadata_sample and vcf_sample columns; requires --metadata')
    args = parser.parse_args(argv)
    if args.id_map and not args.metadata:
        parser.error('--id-map requires --metadata')
    output = Path(args.out)
    if output.exists():
        parser.error(f'Refusing to overwrite {output}')
    try:
        samples, _, _ = vcf_header(args.vcf)
        metadata = read_metadata_groups(args.metadata) if args.metadata else {}
        groups = _resolve_metadata_groups(metadata, samples, args.id_map)
    except (OSError, ValueError, csv.Error) as error:
        parser.error(str(error))
    write_table(output, [{'sample': sample, 'group': groups.get(sample, 'UNRESOLVED')} for sample in samples],
                ['sample', 'group'])
    if args.metadata:
        print(f'Wrote {len(samples)} exact VCF sample IDs to {output}. '
              f'Matched: {len(groups)}; unresolved: {len(samples) - len(groups)}; '
              f'unused metadata: {len(metadata) - len(groups)}.')
        if len(groups) < len(samples):
            print('Replace every UNRESOLVED label using verified metadata.')
    else:
        print(f'Wrote {len(samples)} exact sample IDs to {output}. '
              'Replace every UNRESOLVED label using verified metadata.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
