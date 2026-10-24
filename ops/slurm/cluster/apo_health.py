#!/usr/bin/env python3
"""Print one line per unhealthy apo replicate. Silence means all is well.

Checks the live state CSV each job writes: step must advance between polls,
and energy, temperature and density must stay physical. Also greps the job
log for crash signatures. Previous step counts are kept in a small state file
so a stalled run is distinguishable from a slow one.
"""
import json, math, re, subprocess, sys, time
from pathlib import Path

ROOT = Path("/scratch/users/rsatija/nnrti-mechanisms-git")
APO = ROOT / "results/md_runs/apo"
POSE = ROOT / "results/md_runs/pose_test"
LOGS = ROOT / "logs/apo_md"
POSE_LOGS = ROOT / "logs/pose_test"
STATE = Path.home() / ".apo_health_state.json"

prev = json.loads(STATE.read_text()) if STATE.exists() else {}
cur, alerts = {}, []

# Which replicates Slurm thinks are running: apo_<label>_<rep>
try:
    rows = subprocess.run(["squeue", "-u", "rsatija", "-h", "-t", "R", "-o", "%j|%M"],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          universal_newlines=True, timeout=60).stdout.strip().splitlines()
except Exception as e:
    print(f"ALERT squeue unavailable: {e}")
    sys.exit(0)

def elapsed_s(t):
    """squeue %M: M:SS, H:MM:SS or D-HH:MM:SS."""
    d = 0
    if "-" in t:
        ds, t = t.split("-", 1)
        d = int(ds)
    parts = [int(x) for x in t.split(":")]
    while len(parts) < 3:
        parts.insert(0, 0)
    return d * 86400 + parts[0] * 3600 + parts[1] * 60 + parts[2]


running = {}
for row in rows:
    n, _, el = row.partition("|")
    m = re.match(r"apo_(.+)_(\d+)$", n)
    if m:
        running[(m.group(1), int(m.group(2)))] = (n, elapsed_s(el))
    mp = re.match(r"pt_(pose_.+)$", n)          # pose-independence test jobs
    if mp:
        running[(mp.group(1), 1)] = (n, elapsed_s(el))

for (label, rep), (jobname, age) in sorted(running.items()):
    d = (POSE if label.startswith("pose_") else APO) / label / f"rep_{rep:02d}"
    csvs = list(d.glob("*_md_state.csv"))
    if not csvs:
        if age > 1800:
            alerts.append(f"ALERT {jobname}: RUNNING {age//60} min with no state CSV")
        continue
    try:
        last = subprocess.run(["tail", "-1", str(csvs[0])],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, timeout=30).stdout.strip()
        step, pe, temp, vol, dens, speed = [float(x) for x in last.split(",")]
    except Exception as e:
        alerts.append(f"ALERT {jobname}: cannot parse state CSV ({e})")
        continue

    key = f"{label}_{rep}"
    cur[key] = step
    cur[key + "_t"] = time.time()
    vkey = key + "_vol"
    if vkey in prev and prev[vkey] > 0:
        swing = abs(vol - prev[vkey]) / prev[vkey]
        if swing > 0.03:
            alerts.append(f"ALERT {jobname}: box volume swung {swing*100:.1f}% "
                          f"({prev[vkey]:.0f} -> {vol:.0f} nm^3), barostat unstable")
    cur[vkey] = vol
    if any(math.isnan(v) or math.isinf(v) for v in (pe, temp, dens)):
        alerts.append(f"ALERT {jobname}: NaN/Inf in state  PE={pe} T={temp} rho={dens}")
        continue
    # Only a stall if the step has not moved AND enough wall time has elapsed.
    # Without the clock, re-arming the monitor right after a manual run compares
    # two samples seconds apart and every job looks stalled.
    # A *decrease* means the replicate restarted from scratch (e.g. its
    # checkpoint was moved aside) -- that is not a stall, so just rebaseline.
    # Only an unchanged counter, for long enough, is a stall.
    if key in prev and step == prev[key]:
        age = time.time() - prev.get(key + "_t", 0)
        if age > 240:
            alerts.append(f"ALERT {jobname}: STALLED, step {int(step)} unchanged "
                          f"for {age/60:.0f} min")
    if not (295.0 <= temp <= 307.0):
        alerts.append(f"ALERT {jobname}: temperature {temp:.1f} K out of range")
    if not (0.99 <= dens <= 1.04):
        alerts.append(f"ALERT {jobname}: density {dens:.3f} g/mL out of range")
    if pe > -1.0e6:
        alerts.append(f"ALERT {jobname}: potential energy {pe:.0f} kJ/mol implausibly high")
    if speed < 5.0:
        alerts.append(f"ALERT {jobname}: speed {speed:.1f} ns/day, abnormally slow")

# Crash signatures in logs touched recently
try:
    hits = subprocess.run(
        ["bash", "-lc",
         f"find {LOGS} {POSE_LOGS} -name '*.log' -mmin -30 -newermt '-30 minutes' 2>/dev/null | "
         f"xargs -r grep -lE 'Traceback|Segmentation fault|Particle coordinate is nan|"
         f"Energy is NaN|CUDA error|out of memory' 2>/dev/null | tail -8"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, timeout=90).stdout.split()
    reported = set(prev.get("_logs_seen", []))
    for h in hits:
        if h in reported:
            continue
        reported.add(h)
        sig = subprocess.run(["bash", "-lc",
              f"grep -hoE 'Traceback|Segmentation fault|Particle coordinate is nan|"
              f"Energy is NaN|CUDA error|out of memory' {h} | sort -u | tr '\\n' ' '"],
              stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, timeout=30).stdout.strip()
        alerts.append(f"ALERT log {Path(h).name}: {sig}")
except Exception:
    reported = set(prev.get("_logs_seen", []))

cur["_logs_seen"] = sorted(reported)
STATE.write_text(json.dumps({**prev, **cur}))
for a in alerts:
    print(a)
