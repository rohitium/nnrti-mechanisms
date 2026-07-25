#!/usr/bin/env bash
# Open a persistent SSH ControlMaster to Sherlock (one Duo auth, ~4h persist).
#
# Requires repo-root .env with:
#   SHERLOCK_USERNAME=...
#   SHERLOCK_PASSWORD=...
#
# Usage:
#   bash scripts/sherlock/connect.sh
#   bash scripts/sherlock/connect.sh check    # connect if needed, then module probe
#
# Remote commands after connect:
#   bash scripts/sherlock/remote.sh 'module avail gromacs'

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
ENV_FILE="${ROOT}/.env"
SSH_CTL="${TMPDIR:-/tmp}/nnrti_sherlock_ctl.sock"

if [[ ! -f "${ENV_FILE}" ]]; then
  echo "Missing ${ENV_FILE} (need SHERLOCK_USERNAME and SHERLOCK_PASSWORD)" >&2
  exit 1
fi

# shellcheck disable=SC1090
set -a
source "${ENV_FILE}"
set +a

USER="${SHERLOCK_USERNAME:-${SHERLOCK_USER:-}}"
PASS="${SHERLOCK_PASSWORD:-}"
if [[ -z "${USER}" || -z "${PASS}" ]]; then
  echo "Set SHERLOCK_USERNAME and SHERLOCK_PASSWORD in ${ENV_FILE}" >&2
  exit 1
fi

REMOTE_HOST="${USER}@login.sherlock.stanford.edu"
export SSH_CTL USER PASS REMOTE_HOST

_is_master_up() {
  ssh -S "${SSH_CTL}" -O check "${REMOTE_HOST}" 2>/dev/null
}

_open_master() {
  if _is_master_up; then
    echo "[ssh] ControlMaster already active: ${SSH_CTL}"
    return 0
  fi

  echo "[ssh] Opening ControlMaster to ${REMOTE_HOST}"
  echo "[ssh] Approve Duo Push when prompted on your phone…"

  if ! command -v expect >/dev/null 2>&1; then
    echo "expect not found; falling back to interactive ssh (enter password + Duo manually)" >&2
    ssh -M -S "${SSH_CTL}" -fN \
      -o ControlPersist=4h \
      -o ServerAliveInterval=60 \
      -o ServerAliveCountMax=3 \
      "${REMOTE_HOST}"
    return 0
  fi

  expect <<'EOF'
set timeout 180
spawn ssh -M -S $env(SSH_CTL) -fN \
  -o ControlPersist=4h \
  -o ServerAliveInterval=60 \
  -o ServerAliveCountMax=3 \
  $env(REMOTE_HOST)

expect {
  -re "(?i)password:" {
    send "$env(PASS)\r"
    exp_continue
  }
  -re "(?i)passcode or option" {
    send "1\r"
    exp_continue
  }
  -re "(?i)duo push" {
    exp_continue
  }
  -re "(?i)permission denied" {
    puts stderr "\n[ssh] Authentication failed"
    exit 1
  }
  timeout {
    puts stderr "\n[ssh] Timed out waiting for auth (approve Duo Push?)"
    exit 1
  }
  eof {}
}
EOF

  sleep 1
  if ! _is_master_up; then
    echo "[ssh] ControlMaster not running after auth attempt" >&2
    exit 1
  fi
  echo "[ssh] ControlMaster ready: ${SSH_CTL}"
}

_probe_modules() {
  ssh -S "${SSH_CTL}" "${REMOTE_HOST}" bash -s <<'REMOTE'
set +u
echo "=== hostname ==="
hostname
echo "=== gromacs (spider) ==="
module spider gromacs 2>&1 || true
echo "=== gromacs (avail) ==="
module avail gromacs 2>&1 || true
echo "=== schrodinger (spider) ==="
module spider schrodinger 2>&1 || true
module spider Schrodinger 2>&1 || true
echo "=== openmm (avail) ==="
module avail py-openmm 2>&1 || true
module spider openmm 2>&1 || true
echo "=== cuda (head) ==="
module avail cuda 2>&1 | head -30 || true
echo "=== python (head) ==="
module avail python 2>&1 | head -20 || true
REMOTE
}

_open_master

case "${1:-}" in
  check|modules|probe)
    _probe_modules
    ;;
  "")
    echo "Use: bash scripts/sherlock/remote.sh '<command>'"
    ;;
  *)
    ssh -S "${SSH_CTL}" "${REMOTE_HOST}" "$@"
    ;;
esac
