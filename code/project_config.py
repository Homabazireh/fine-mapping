"""Load validated project settings using only the Python standard library."""
from __future__ import annotations

import argparse
from decimal import Decimal, InvalidOperation
import json
import math
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CONFIG = PROJECT_ROOT / 'config' / 'settings.json'
DEFAULT_PATHS = {
    'vcf': 'data/raw/full_LR-WGS.vcf.gz',
    'gwas': 'data/raw/EUR_diagram_t2d_1e-8.tsv',
    'metadata': 'data/raw/participant_categories.tsv',
    'samples': 'data/intermediate/samples.tsv',
    'out': 'results/cohort_categories',
    'demo_data': 'data/intermediate/demo',
    'demo_results': 'results/demo',
}
DEFAULT_ANALYSIS = {
    'p_threshold': '1e-8',
    'lead_gap': 150000,
    'top_loci': 20,
    'flank': 50000,
    'max_markers': 8,
    'min_dp': 10.0,
    'min_gq': 20.0,
    'min_call_rate': 0.90,
    'min_mac': 2,
    'max_indel_length': 49,
    'allow_filter_dot': False,
    'allow_subset': False,
    'expected_per_group': None,
    'test_associations': False,
    'min_test_samples': 10,
    'min_resolved_fraction': 0.80,
    'min_carriers': 3,
    'plot': False,
    'max_plots': 10,
}


def resolve_project_path(path):
    """Resolve a settings path relative to the project, preserving raw symlinks."""
    path = Path(path).expanduser()
    return path if path.is_absolute() else PROJECT_ROOT / path


def _unique_object(pairs):
    obj = {}
    for key, value in pairs:
        if key in obj:
            raise ValueError(f'Duplicate JSON setting: {key}')
        obj[key] = value
    return obj


def _validate_analysis(settings):
    for key, value in settings.items():
        default = DEFAULT_ANALYSIS[key]
        if key == 'p_threshold':
            if isinstance(value, bool) or not isinstance(value, (str, int, float)):
                raise ValueError('analysis.p_threshold must be a decimal string or number')
            try:
                threshold = Decimal(str(value))
            except InvalidOperation as error:
                raise ValueError('analysis.p_threshold must be a valid decimal') from error
            if not threshold.is_finite() or not 0 <= threshold <= 1:
                raise ValueError('analysis.p_threshold must be within [0,1]')
            settings[key] = str(value)
        elif isinstance(default, bool):
            if not isinstance(value, bool):
                raise ValueError(f'analysis.{key} must be true or false')
        elif key == 'expected_per_group':
            if value is not None and (type(value) is not int or value < 1):
                raise ValueError('analysis.expected_per_group must be null or a positive integer')
        elif isinstance(default, int):
            if type(value) is not int:
                raise ValueError(f'analysis.{key} must be an integer')
            minimum = 2 if key == 'max_markers' else (1 if key in {
                'max_indel_length', 'min_test_samples', 'min_carriers', 'max_plots'} else 0)
            if value < minimum or (key == 'max_indel_length' and value > 49):
                raise ValueError(f'analysis.{key} is outside its allowed range')
        else:
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError(f'analysis.{key} must be a finite number')
            if value < 0 or (key in {'min_call_rate', 'min_resolved_fraction'} and value > 1):
                raise ValueError(f'analysis.{key} is outside its allowed range')
            if key == 'min_call_rate' and value == 0:
                raise ValueError('analysis.min_call_rate must be greater than zero')


def load_config(config_path=None):
    """Read settings, merge built-in defaults, and return absolute path strings.

    Custom files may override a subset of settings. Both their relative path
    values and a relative config filename use PROJECT_ROOT as their base.
    Scientific attestations and phenotype mappings are intentionally not settings.
    """
    path = resolve_project_path(config_path) if config_path is not None else DEFAULT_CONFIG
    with path.open(encoding='utf-8-sig') as handle:
        data = json.load(handle, object_pairs_hook=_unique_object)
    if not isinstance(data, dict):
        raise ValueError(f'{path}: settings must be a JSON object')
    unexpected = set(data) - {'paths', 'analysis'}
    if unexpected:
        raise ValueError(f'{path}: unknown settings sections: {", ".join(sorted(unexpected))}')
    merged = {}
    for section, defaults in [('paths', DEFAULT_PATHS), ('analysis', DEFAULT_ANALYSIS)]:
        updates = data.get(section, {})
        if not isinstance(updates, dict):
            raise ValueError(f'{path}: {section} must be a JSON object')
        unexpected = set(updates) - set(defaults)
        if unexpected:
            raise ValueError(f'{path}: unknown {section} settings: {", ".join(sorted(unexpected))}')
        merged[section] = {**defaults, **updates}
    for key, value in merged['paths'].items():
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f'{path}: paths.{key} must be a nonempty path string')
        merged['paths'][key] = str(resolve_project_path(value))
    _validate_analysis(merged['analysis'])
    return merged


class ConfigArgumentParser(argparse.ArgumentParser):
    """Apply project defaults before parsing explicit command-line overrides."""

    def __init__(self, *args, path_defaults=None, analysis_defaults=False, **kwargs):
        kwargs.setdefault('allow_abbrev', False)
        super().__init__(*args, **kwargs)
        self.path_defaults = dict(path_defaults or {})
        self.analysis_defaults = analysis_defaults
        self.add_argument('--config', default=str(DEFAULT_CONFIG),
                          help='Settings JSON; relative config paths use the project root')

    def parse_known_args(self, args=None, namespace=None):
        selector = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
        selector.add_argument('--config', default=str(DEFAULT_CONFIG))
        selected, _ = selector.parse_known_args(args)
        try:
            settings = load_config(selected.config)
        except (OSError, ValueError) as error:
            self.error(str(error))
        defaults = {argument: settings['paths'][setting]
                    for argument, setting in self.path_defaults.items()}
        if self.analysis_defaults:
            defaults.update(settings['analysis'])
        self.set_defaults(**defaults)
        parsed, remaining = super().parse_known_args(args, namespace)
        parsed.config = str(resolve_project_path(parsed.config))
        return parsed, remaining
