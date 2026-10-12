#!/usr/bin/env python3
"""Run all project tests from any working directory."""
from pathlib import Path
import sys
import unittest


if __name__ == "__main__":
    root = Path(__file__).resolve().parent.parent
    sys.path.insert(0, str(root / "code"))
    suite = unittest.defaultTestLoader.discover(
        str(root / "tests"), pattern="test_*.py"
    )
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    raise SystemExit(0 if result.wasSuccessful() else 1)
