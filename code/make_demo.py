#!/usr/bin/env python3
"""Create a small, entirely synthetic test cohort. NOT participant data or results."""
from pathlib import Path
from project_config import ConfigArgumentParser


def make_demo(folder):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    samples = [f'{prefix}{i}' for prefix in ['D', 'N', 'C'] for i in range(1, 5)]
    groups = ['T2D'] * 4 + ['Prediabetes'] * 4 + ['NORM'] * 4
    with open(folder / 'samples.tsv', 'w') as f:
        f.write('sample\tgroup\n')
        for sample, group in zip(samples, groups):
            f.write(f'{sample}\t{group}\n')
    (folder / 'gwas.tsv').write_text(
        'Chromsome\tPosition\tEffectAllele\tNonEffectAllele\tBeta\tSE\tEAF\tPval\tNcases\tNcontrols\tNeff\tliftover_strand\n'
        '1\t10000\tG\tA\t0.2\t0.03\t0.3\t1e-12\t1000\t1000\t2000\t+\n'
        '1\t9995\tA\tAT\t0.1\t0.03\t0.2\t1e-9\t1000\t1000\t2000\t+\n'
        '1\t200010\tT\tC\t0.2\t0.04\t0.4\t1e-10\t1000\t1000\t2000\t+\n'
        '1\t200018\tA\tG\t0.1\t0.03\t0.2\t1e-8\t1000\t1000\t2000\t+\n')
    header = ('##fileformat=VCFv4.2\n##reference=GRCh38\n##contig=<ID=chr1,length=1000000>\n'
              '##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate">\n'
              '##FILTER=<ID=PASS,Description="All filters passed">\n'
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
              '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n'
              '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n'
              '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n'
              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(samples) + '\n')
    haplotypes = [(0, 0, 0, 0), (1, 1, 1, 1), (1, 0, 1, 0)]
    pairs = [(1, 1), (1, 0), (1, 2), (1, 0), (1, 0), (1, 2), (0, 0), (0, 2),
             (0, 0), (0, 0), (0, 2), (0, 0)]
    records = []
    for j, (pos, ref, alt) in enumerate([(9990, 'C', 'T'), (9995, 'AT', 'A'), (10000, 'A', 'G'), (10005, 'G', 'C')]):
        cells = []
        for i, sample in enumerate(samples):
            h1, h2 = pairs[i]
            a, b = haplotypes[h1][j], haplotypes[h2][j]
            ps = 1000 + i
            if sample == 'D4' and pos == 10005:
                ps += 100  # A real PS boundary in the synthetic fixture: do not stitch.
            gt = f'{a}|{b}:{25}:{40}:{ps}'
            if sample == 'N4' and pos == 9995:
                gt = './.:.:.:.'  # Must not become 0/0.
            cells.append(gt)
        records.append(f'chr1\t{pos}\t.\t{ref}\t{alt}\t60\tPASS\t.\tGT:DP:GQ:PS\t' + '\t'.join(cells))
    # Deliberately out-of-scope records test exclusion; no repeat analysis is performed.
    for pos, ref, alt, info in [(10030, 'A', '<DEL>', 'END=10100'), (10040, 'A', 'A' * 55, '.')]:
        records.append(f'chr1\t{pos}\t.\t{ref}\t{alt}\t60\tPASS\t{info}\tGT:DP:GQ:PS\t' + '\t'.join(['0|1:25:40:1000'] * 12))
    for j, (pos, ref, alt) in enumerate([(200005, 'T', 'TA'), (200010, 'C', 'T'), (200018, 'G', 'A')]):
        cells = []
        for i in range(12):
            a, b = (i % 2), ((i // 2) % 2)
            cells.append(f'{a}|{b}:25:40:{2000 + i}')
        records.append(f'chr1\t{pos}\t.\t{ref}\t{alt}\t60\tPASS\t.\tGT:DP:GQ:PS\t' + '\t'.join(cells))
    (folder / 'cohort.vcf').write_text(header + '\n'.join(records) + '\n')
    (folder / 'README.txt').write_text('All records, participants, phenotypes and GWAS statistics here are SYNTHETIC software-test fixtures.\n')
    return folder


def main(argv=None):
    parser = ConfigArgumentParser(description=__doc__, path_defaults={'out': 'demo_data'})
    parser.add_argument('--out', type=Path,
                        help='Directory for synthetic fixtures (default: configured demo_data path)')
    args = parser.parse_args(argv)
    folder = make_demo(args.out)
    print(f'Wrote entirely synthetic test fixtures to {folder}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
