#!/usr/bin/env python3
"""Run a synthetic demonstration with its own VCF, GWAS, and sample table."""
import argparse
from pathlib import Path
import tempfile

from compare_haplotypes import main as run_analysis
from make_demo import make_demo
from project_config import DEFAULT_CONFIG, load_config


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('demo_folder', nargs='?', type=Path,
                        help='Override the configured synthetic input directory')
    parser.add_argument('output_folder', nargs='?', type=Path,
                        help='New or empty output directory; default: a new folder under results/demo')
    parser.add_argument('--config', default=str(DEFAULT_CONFIG),
                        help='Project settings JSON; relative config paths are project-relative')
    args = parser.parse_args(argv)
    try:
        settings = load_config(args.config)
        fixtures = make_demo(args.demo_folder or settings['paths']['demo_data'])
        output = args.output_folder
        if output is None:
            parent = Path(settings['paths']['demo_results'])
            parent.mkdir(parents=True, exist_ok=True)
            output = Path(tempfile.mkdtemp(prefix='run_', dir=parent))
    except (ValueError, OSError) as error:
        parser.error(str(error))
    return run_analysis([
        '--config', str(args.config),
        '--vcf', str(fixtures / 'cohort.vcf'),
        '--gwas', str(fixtures / 'gwas.tsv'),
        '--samples', str(fixtures / 'samples.tsv'),
        '--out', str(output),
        '--genome-build', 'GRCh38',
        '--gwas-alleles-forward-confirmed',
        '--expected-per-group', '4',
    ])


if __name__ == '__main__':
    raise SystemExit(main())
