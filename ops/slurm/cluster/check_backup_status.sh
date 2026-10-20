#!/bin/bash
#
# Report the state of the FEP raw-data backup job.
#
# Answers three questions in one shot:
#   1. Was the job submitted?   -- squeue / sacct
#   2. Is it running?           -- STATUS.txt state= line
#   3. Has it finished?         -- STATUS.txt state=DONE plus the manifest
#
# Runs from either side. On the Mac it goes over the existing SSH connection;
# on a Sherlock login node it runs the same checks locally. Read-only: it
# submits nothing and changes nothing.
#
# USAGE
#   bash ops/slurm/cluster/check_backup_status.sh            # today's archive
#   STAMP=2026-10-20 bash ops/slurm/cluster/check_backup_status.sh
#   LINES=40 bash ops/slurm/cluster/check_backup_status.sh    # longer log tail
#
set -euo pipefail

STAMP="${STAMP:-$(date +%Y-%m-%d)}"
LINES="${LINES:-15}"
SSH_HOST="${SSH_HOST:-sherlock}"

read -r -d '' BODY <<BODYEOF || true
DEST=\$GROUP_HOME/nnrti-mechanisms/fep_raw/$STAMP

echo '=== Slurm queue ==='
squeue -u \$USER -o '%.10i %.14j %.9T %.10M %.10l %R' 2>/dev/null \
    | grep -E 'JOBID|fep_backup|fep_verify' || echo '  (no backup or verify job in the queue)'

echo
echo '=== Recent jobs (sacct, last 2 days) ==='
sacct -u \$USER --starttime now-2days \
    --format=JobID%14,JobName%16,State%12,Elapsed%10 2>/dev/null \
    | grep -E 'JobID|---|fep_backup|fep_verify' | head -14 || echo '  (sacct unavailable)'

echo
echo "=== STATUS (\$DEST) ==="
if [ -f "\$DEST/STATUS.txt" ]; then cat "\$DEST/STATUS.txt"; else echo '  (no STATUS.txt yet)'; fi

echo
echo "=== PROGRESS.log (last $LINES lines) ==="
if [ -f "\$DEST/PROGRESS.log" ]; then tail -n $LINES "\$DEST/PROGRESS.log"; else echo '  (no PROGRESS.log yet)'; fi

echo
echo '=== Archives ==='
if [ -d "\$DEST" ]; then
    n=\$(ls -1 "\$DEST"/*.tar.gz 2>/dev/null | wc -l)
    echo "  \$n archives, \$(du -sh "\$DEST" 2>/dev/null | cut -f1) total"
    for f in "\$DEST"/*.part "\$DEST"/*.CORRUPT "\$DEST"/*.unverified; do
        [ -e "\$f" ] && echo "  INCOMPLETE: \$f"
    done
    echo "  manifest rows: \$(( \$(wc -l < "\$DEST/BACKUP_MANIFEST.csv" 2>/dev/null || echo 1) - 1 ))"
else
    echo '  (destination not created yet)'
fi
BODYEOF

if [ -n "${GROUP_HOME:-}" ] && [ -d "${GROUP_HOME:-/nonexistent}" ]; then
    # Already on Sherlock.
    bash -c "$BODY"
else
    ssh -o BatchMode=yes "$SSH_HOST" "$BODY"
fi
