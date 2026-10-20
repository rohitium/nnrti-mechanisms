#!/bin/bash
#
# Report the state of the FEP raw-data backup job on Sherlock.
#
# Answers three questions in one shot:
#   1. Was the job submitted?   -- squeue / sacct
#   2. Is it running?           -- STATUS.txt state= line
#   3. Has it finished?         -- STATUS.txt state=DONE plus the manifest
#
# Runs from the Mac over the existing SSH connection. Read-only: it submits
# nothing and changes nothing on the cluster.
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

ssh -o BatchMode=yes "$SSH_HOST" "
DEST=\$GROUP_HOME/nnrti-mechanisms/fep_raw/$STAMP

echo '=== Slurm queue (fep_backup jobs) ==='
squeue -u \$USER -o '%.10i %.14j %.9T %.10M %.10l %R' 2>/dev/null | grep -E 'JOBID|fep_backup' || echo '  (no fep_backup job in the queue)'

echo
echo '=== Recent fep_backup jobs (sacct, last 2 days) ==='
sacct -u \$USER --starttime now-2days --name=fep_backup_t3,fep_backup_t2,fep_backup_t1 \
    --format=JobID%12,JobName%16,State%12,Elapsed%10,MaxRSS%10 2>/dev/null | head -12 \
    || echo '  (sacct unavailable)'

echo
echo \"=== STATUS (\$DEST) ===\"
if [ -f \"\$DEST/STATUS.txt\" ]; then cat \"\$DEST/STATUS.txt\"; else echo '  (no STATUS.txt yet -- job has not started writing)'; fi

echo
echo \"=== PROGRESS.log (last $LINES lines) ===\"
if [ -f \"\$DEST/PROGRESS.log\" ]; then tail -n $LINES \"\$DEST/PROGRESS.log\"; else echo '  (no PROGRESS.log yet)'; fi

echo
echo '=== Archives written so far ==='
if [ -d \"\$DEST\" ]; then
    ls -lh \"\$DEST\"/*.tar.gz 2>/dev/null | awk '{print \"  \", \$9, \$5}' || echo '  (none yet)'
    echo -n '  total: '; du -sh \"\$DEST\" 2>/dev/null | cut -f1
else
    echo '  (destination not created yet)'
fi
"
