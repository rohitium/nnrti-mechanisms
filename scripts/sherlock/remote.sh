#!/usr/bin/env bash
# Run a command on Sherlock via the ControlMaster from connect.sh
#
# Usage:
#   bash scripts/sherlock/remote.sh 'module avail gromacs'
#   bash scripts/sherlock/remote.sh hostname

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
ENV_FILE="${ROOT}/.env"
SSH_CTL="${TMPDIR:-/tmp}/nnrti_sherlock_ctl.sock"

if [[ ! -f "${ENV_FILE}" ]]; then
  echo "Missing ${ENV_FILE}" >&2
  exit 1
fi

set -a
# shellcheck disable=SC1090
source "${ENV_FILE}"
set +a

USER="${SHERLOCK_USERNAME:-${SHERLOCK_USER:-}}"
REMOTE_HOST="${USER}@login.sherlock.stanford.edu"

if [[ $# -lt 1 ]]; then
  echo "Usage: bash scripts/sherlock/remote.sh '<remote command>'" >&2
  exit 1
fi

if ! ssh -S "${SSH_CTL}" -O check "${REMOTE_HOST}" 2>/dev/null; then
  echo "No ControlMaster — run: bash scripts/sherlock/connect.sh" >&2
  exit 1
fi

ssh -S "${SSH_CTL}" "${REMOTE_HOST}" "$@"
