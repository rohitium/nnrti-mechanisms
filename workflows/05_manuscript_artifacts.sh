#!/usr/bin/env bash
# Export exact submitted artwork/tables and rebuild current numerical artifacts.
set -euo pipefail
cd "$(dirname "$0")/.."
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
exec "${PYTHON:-python}" -m nnrti.cli.manuscript --export --rebuild all "$@"
