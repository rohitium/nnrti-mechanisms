#!/usr/bin/env python3
"""Repair classical-MD run metadata to match ground truth. DRY-RUN by default.

WARNING: ``state.csv`` is *not* always ground truth for the holo 100 ns panel.
Many analysis DCDs match the completed 100 ns fingerprint (180 frames / 34 MB
or 360 frames / 68 MB) while JSON+CSV still say 50–70 ns — leftover mid-slice
``status=ok`` dumps. Do not blindly copy those step counts into the JSON.
Use ``nnrti.analysis.md_timing.infer_production_ns`` (DCD fingerprint) for
analysis time axes. This script only reconciles JSON↔CSV agreement.

Nothing is written unless ``--apply`` is passed. All actions are idempotent.

Actions:
  1. Archive superseded / stub metadata to ``results/md_runs/_archive/`` and log
     it in ``manifests/md_archive_log.csv`` (never delete):
       * ``*_quick.json`` short-test metadata (wt/rep_01, Y188L/rep_01)
       * the aborted 4000-step stub ``apo/y188l/rep_01`` (light files only)
  2. Fix stale stored paths in apo JSONs (results/apo_md_runs -> results/md_runs/apo),
     verifying the corrected path exists on disk.
  3. Correct ``md_production_steps_completed`` (and ``md_total_steps``) where the
     JSON disagrees with ``state.csv`` (e.g. apo/k103n_m230l).
  4. Label ``results/md_runs/wt_sherlock_rerun`` as a diagnostic, non-canonical set.

    python3 ops/maintenance/md/repair_md_metadata.py            # dry-run: print every change
    python3 ops/maintenance/md/repair_md_metadata.py --apply    # perform the changes
"""
from __future__ import annotations

import argparse
import csv
import datetime as _dt
import hashlib
import json
import shutil
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from nnrti.paths import MD_RUNS, APO_MD_RUNS, MANIFESTS, REPO_ROOT as ROOT, rel  # noqa: E402

DT_FS = 2.0
STALE = "apo_md_runs"
CANON = "md_runs/apo"
APO_TARGET_STEPS = 50_000_000  # canonical 100 ns production target for apo runs
ARCHIVE_DIR = MD_RUNS / "_archive"
ARCHIVE_LOG = MANIFESTS / "md_archive_log.csv"


def _sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _last_state_step(rundir: Path) -> int | None:
    for csvf in sorted(rundir.glob("*_md_state.csv")):
        last = None
        with open(csvf) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                last = line
        if last:
            try:
                return int(float(last.split(",")[0]))
            except ValueError:
                pass
    return None


class Planner:
    def __init__(self, apply: bool):
        self.apply = apply
        self.log_rows: list[dict] = []
        self.archived: set[Path] = set()
        self.n = 0

    def _log_archive(self, original: Path, archived: Path, reason: str):
        self.log_rows.append({
            "date": _dt.date.today().isoformat(),
            "reason": reason,
            "original_path": rel(original),
            "archived_path": rel(archived),
            "sha256": _sha256(archived if archived.exists() else original),
            "bytes": (archived if archived.exists() else original).stat().st_size,
        })

    def archive(self, f: Path, reason: str):
        dest = ARCHIVE_DIR / f.relative_to(MD_RUNS)
        self.n += 1
        self.archived.add(f.resolve())
        print(f"[archive] {rel(f)}\n          -> {rel(dest)}   ({reason})")
        if self.apply:
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(f), str(dest))
            self._log_archive(f, dest, reason)

    def fix_paths(self, jf: Path):
        d = json.loads(jf.read_text())
        changes = {}
        for k, v in d.items():
            if isinstance(v, str) and STALE in v:
                nv = v.replace(f"results/{STALE}", "results/" + CANON)
                exists = (ROOT / nv).exists()
                changes[k] = (v, nv, exists)
        if not changes:
            return
        self.n += 1
        print(f"[fix_path] {rel(jf)}")
        for k, (old, new, exists) in changes.items():
            flag = "" if exists else "   <-- WARNING: target missing"
            print(f"           {k}: {old}\n              -> {new}{flag}")
        if self.apply:
            for k, (_old, new, _e) in changes.items():
                d[k] = new
            jf.write_text(json.dumps(d, indent=2) + "\n")

    def fix_steps(self, jf: Path, state_step: int):
        d = json.loads(jf.read_text())
        old_completed = d.get("md_production_steps_completed")
        if old_completed == state_step:
            return
        old_total = d.get("md_total_steps")
        heating = (old_total - old_completed) if (old_total and old_completed is not None) else 12500
        new_total = heating + state_step
        self.n += 1
        print(f"[fix_steps] {rel(jf.parent)}")
        print(f"            md_production_steps_completed: {old_completed} -> {state_step}"
              f"   ({state_step * DT_FS / 1e6:.1f} ns from state.csv)")
        print(f"            md_total_steps:                {old_total} -> {new_total}")
        if self.apply:
            d["md_production_steps_completed"] = state_step
            d["md_total_steps"] = new_total
            d["production_ns_completed"] = round(state_step * DT_FS / 1e6, 3)
            d["metadata_repaired"] = _dt.date.today().isoformat()
            jf.write_text(json.dumps(d, indent=2) + "\n")

    def fix_apo_target(self, jf: Path):
        """Normalize the apo production *target* to 100 ns.

        The original apo JSONs recorded a 5M-step (10 ns) target, but the runs
        were intended for 100 ns (config default; matches the manual Sherlock
        edits). Left uncorrected, a run that reached 18M steps reads as
        'completed 18M > target 5M', which is nonsensical.
        """
        d = json.loads(jf.read_text())
        old = d.get("md_production_steps")
        if old is None or old == APO_TARGET_STEPS:
            return
        self.n += 1
        print(f"[apo_target] {rel(jf.parent)}: md_production_steps {old} -> {APO_TARGET_STEPS} (100 ns intent)")
        if self.apply:
            d["md_production_steps"] = APO_TARGET_STEPS
            jf.write_text(json.dumps(d, indent=2) + "\n")

    def label_diagnostic(self, rundir: Path, text: str):
        note = rundir / "DIAGNOSTIC.md"
        if note.exists():
            return
        self.n += 1
        print(f"[label] {rel(note)} (mark non-canonical diagnostic set)")
        if self.apply:
            note.write_text(text)

    def flush_log(self):
        if not (self.apply and self.log_rows):
            return
        ARCHIVE_LOG.parent.mkdir(parents=True, exist_ok=True)
        exists = ARCHIVE_LOG.exists()
        with open(ARCHIVE_LOG, "a", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["date", "reason", "original_path", "archived_path", "sha256", "bytes"])
            if not exists:
                w.writeheader()
            for row in self.log_rows:
                w.writerow(row)
        print(f"\n[log] appended {len(self.log_rows)} row(s) to {rel(ARCHIVE_LOG)}")


def build_plan(p: Planner):
    # 1. archive superseded quick-test metadata
    for q in sorted(MD_RUNS.glob("*/rep_*/*_quick.json")):
        p.archive(q, "superseded quick-test metadata")

    # 1b. archive the aborted apo/y188l/rep_01 stub (light metadata only)
    stub = APO_MD_RUNS / "y188l" / "rep_01"
    if stub.exists() and (_last_state_step(stub) or 0) < 50_000:
        for f in sorted(list(stub.glob("*.json")) + list(stub.glob("*_md_state.csv"))):
            p.archive(f, "aborted 4000-step (8 ps) test stub")

    # 2. fix stale apo paths
    for jf in sorted(APO_MD_RUNS.rglob("*.json")):
        if not jf.exists() or jf.resolve() in p.archived:
            continue
        try:
            txt = jf.read_text()
        except Exception:
            continue
        if STALE in txt:
            p.fix_paths(jf)

    # 3. fix step counts where JSON disagrees with state.csv (holo + apo)
    rep_dirs = [d for d in MD_RUNS.glob("*/rep_*") if d.is_dir() and d.parent.name != "apo"]
    rep_dirs += [d for d in APO_MD_RUNS.glob("*/rep_*") if d.is_dir()]
    for rd in sorted(rep_dirs):
        state = _last_state_step(rd)
        if state is None:
            continue
        for jf in sorted(rd.glob("*.json")):
            if jf.resolve() in p.archived:
                continue
            try:
                d = json.loads(jf.read_text())
            except Exception:
                continue
            if "md_production_steps" not in d:
                continue
            if d.get("md_production_steps_completed") not in (None, state):
                p.fix_steps(jf, state)

    # 3b. normalize apo production target to 100 ns (fixes completed > target)
    for jf in sorted(APO_MD_RUNS.rglob("*.json")):
        if not jf.exists() or jf.resolve() in p.archived:
            continue
        try:
            d = json.loads(jf.read_text())
        except Exception:
            continue
        if "md_production_steps" in d:
            p.fix_apo_target(jf)

    # 4. label the diagnostic WT rerun
    rerun = MD_RUNS / "wt_sherlock_rerun"
    if rerun.exists():
        lengths = {rd.name: _last_state_step(rd) for rd in sorted(rerun.glob("rep_*"))}
        ns = {k: (v * DT_FS / 1e6 if v else None) for k, v in lengths.items()}
        text = (
            "# wt_sherlock_rerun — DIAGNOSTIC, non-canonical\n\n"
            "Purpose: repeat WT holo MD to assess whether the wide spread in MM/GBSA\n"
            "energies was statistical noise vs. a single simulation-artifact outlier.\n\n"
            "**Not part of the canonical panel.** Canonical WT = `results/md_runs/wt/`.\n\n"
            f"Status: incomplete — per-replicate length from state.csv: {ns}\n"
            "No run JSON was written (a limitation this metadata cleanup documents).\n"
        )
        p.label_diagnostic(rerun, text)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Repair MD run metadata from ground truth (dry-run by default).")
    ap.add_argument("--apply", action="store_true", help="Perform changes (default: dry-run).")
    args = ap.parse_args(argv)

    mode = "APPLY" if args.apply else "DRY-RUN (no files written)"
    print(f"=== MD metadata repair — {mode} ===\n")
    p = Planner(apply=args.apply)
    build_plan(p)
    p.flush_log()
    print(f"\n{'Applied' if args.apply else 'Planned'} {p.n} change(s).")
    if not args.apply and p.n:
        print("Re-run with --apply to perform them.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
