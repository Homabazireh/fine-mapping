"""Readable group names for figures and reports, without changing assignments.

The participant table remains the source of group membership. These helpers
only format its category labels; they never combine groups or reassign people.
"""


_DISPLAY_NAMES = {
    'T2D': 'Diabetes',
    'diabetes': 'Diabetes',
    'Diabetes': 'Diabetes',
    'Prediabetes': 'Prediabetes',
    'prediabetes': 'Prediabetes',
    'NORM': 'Control',
    'control': 'Control',
    'Control': 'Control',
}
_GROUP_ORDER = {'Diabetes': 0, 'Prediabetes': 1, 'Control': 2}


def display_group(group):
    """Expand known source labels; leave any unrecognized category unchanged."""
    return _DISPLAY_NAMES.get(group, group)


def ordered_groups(groups):
    """Return original group keys in display order, retaining separate groups."""
    return sorted(groups, key=lambda group: (
        _GROUP_ORDER.get(display_group(group), len(_GROUP_ORDER)), str(group),
    ))
