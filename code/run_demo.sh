#!/usr/bin/env bash
# Synthetic inputs and outputs are kept in data/intermediate and results.
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
exec python3 "$script_dir/run_demo.py" "$@"
