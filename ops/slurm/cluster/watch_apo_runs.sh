#!/bin/bash
# Emit one line per actionable change in the apo queue. Silence = steady state.
SINCE="${SINCE:-$(date -u -v-1d +%Y-%m-%dT%H:%M:%S 2>/dev/null || date -u -d '1 day ago' +%Y-%m-%dT%H:%M:%S)}"
prev_sig=""
# Reported failures persist across re-arms, otherwise every restart of this
# watcher re-reports every failure inside the SINCE window.
SEEN_FILE="${SEEN_FILE:-$HOME/.apo_watch_seen}"
touch "$SEEN_FILE"
while true; do
  q=$(ssh -o BatchMode=yes -o ConnectTimeout=30 sherlock 'squeue -u $USER -h -o "%j|%T|%r"' 2>/dev/null)
  if [ -z "$q" ]; then
    bad=$(ssh -o BatchMode=yes -o ConnectTimeout=30 sherlock "sacct -u \$USER -S $SINCE -X -P -n --format=JobID,JobName,State" 2>/dev/null \
          | grep -E "FAILED|TIMEOUT|NODE_FAIL" | grep -E "apo_|pt_pose")
    [ -n "$bad" ] && echo "QUEUE EMPTY. terminal failures since $SINCE:" && echo "$bad"
    echo "QUEUE EMPTY - all apo jobs finished or died; check and release held work"
    exit 0
  fi
  run=$(printf '%s\n' "$q" | grep -c '|RUNNING|')
  elig=$(printf '%s\n' "$q" | grep '|PENDING|' | grep -vc 'JobHeldUser')
  held=$(printf '%s\n' "$q" | grep -c 'JobHeldUser')

  # Terminal failures we have not reported yet
  bad=$(ssh -o BatchMode=yes -o ConnectTimeout=30 sherlock "sacct -u \$USER -S $SINCE -X -P -n --format=JobID,JobName,State,Elapsed" 2>/dev/null \
        | grep -E "FAILED|TIMEOUT|NODE_FAIL" | grep -E "apo_|pt_pose")
  while IFS= read -r line; do
    [ -z "$line" ] && continue
    id=${line%%|*}
    grep -qx "$id" "$SEEN_FILE" 2>/dev/null && continue
    echo "JOB FAILED: $line"
    echo "$id" >> "$SEEN_FILE"
  done <<< "$bad"

  sig="R=$run ELIG=$elig HELD=$held"
  if [ "$sig" != "$prev_sig" ]; then
    [ -n "$prev_sig" ] && echo "STATE $sig (was ${prev_sig})"
    prev_sig="$sig"
  fi
  # Capacity free and priority queue drained -> time to release held work
  if [ "$elig" -eq 0 ] && [ "$held" -gt 0 ] && [ "$run" -lt 16 ]; then
    echo "CAPACITY FREE: running=$run, no eligible jobs, $held held - release next batch"
  fi
  # Per-replicate health: step advancing, energy/temperature/density physical,
  # logs free of crash signatures. Prints nothing when everything is fine.
  ssh -o BatchMode=yes -o ConnectTimeout=30 sherlock 'python3 ~/apo_health.py' 2>/dev/null

  sleep 300
done
