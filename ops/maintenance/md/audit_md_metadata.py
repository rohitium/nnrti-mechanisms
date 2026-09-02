#!/usr/bin/env python3
"""Audit classical-MD run metadata for consistency — read-only.

For every run under ``results/md_runs`` (holo) and ``results/md_runs/apo`` (apo),
cross-check the run JSON against the ground-truth signals actually on disk:

  * ``state.csv`` last Step      — what OpenMM's StateDataReporter wrote
  * checkpoint ``currentStep``   — what a resume literally starts from (``--check-checkpoints``; needs OpenMM)
  * presence of ``.chk`` + ``assets/*_system.xml`` — is the run resumable?

and flag every inconsistency class:

  * JSON step count disagreeing with ``state.csv`` (silently stale metadata),
  * stale stored paths in the JSON (e.g. old ``results/apo_md_runs/...``),
  * duplicate or stub JSONs in a run directory.

JSON is a claim, never trusted. ``state.csv`` is better but can also be a stale
mid-slice dump. For analysis time axes, prefer the analysis-DCD fingerprint
(see ``nnrti.analysis.md_timing``): 180 fr / 34 MB and 360 fr / 68 MB are the
completed 100 ns products, even when JSON/CSV say 50–70 ns.

    python3 ops/maintenance/md/audit_md_metadata.py                    # table + summary
    python3 ops/maintenance/md/audit_md_metadata.py --check-checkpoints # also read .chk currentStep (needs openmm)
    python3 ops/maintenance/md/audit_md_metadata.py --format csv > md_metadata_audit.csv
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import dataclass, asdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from nnrti.paths import MD_RUNS, APO_MD_RUNS, REPO_ROOT as _ROOT, rel  # noqa: E402

DEFAULT_TIMESTEP_FS = 2.0
TARGET_1US_STEPS = 500_000_000  # 1 us at 2 fs, for "remaining" reporting


@dataclass
class RunAudit:
    run: str
    phase: str
    json_files: int          # how many *.json with a step field live in the dir
    json_steps: int | None   # md_production_steps_completed from the primary JSON
    state_steps: int | None  # last Step in *_md_state.csv (ground truth)
    ckpt_steps: int | None    # currentStep read from the .chk (if --check-checkpoints)
    ns_done: float | None
    json_matches_state: bool | None
    ckpt_matches_state: bool | None
    stored_path_ok: bool | None
    has_chk: bool
    has_system_xml: bool
    status: str | None
    diagnostic: bool = False

    @property
    def has_data(self) -> bool:
        """A run with any actual output (state log or metadata). Empty/reserved
        placeholder dirs have neither and are not inconsistencies."""
        return self.state_steps is not None or self.json_files > 0

    @property
    def clean(self) -> bool:
        if not self.has_data or self.diagnostic:
            return True  # empty/reserved dir, or a labeled non-canonical diagnostic set
        issues = [
            self.json_matches_state is False,
            self.stored_path_ok is False,
            self.json_files != 1,
            self.ckpt_matches_state is False,
        ]
        return not any(issues)


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


def _checkpoint_step(chk: Path) -> int | None:
    """Read OpenMM checkpoint currentStep without running dynamics. Needs openmm."""
    try:
        import openmm  # noqa: F401
        from openmm import app
        import openmm as mm
        import openmm.unit as unit
    except Exception:
        return None
    try:
        # A checkpoint needs a matching System+topology to load into a Context.
        assets = list((chk.parent / "assets").glob("*_system.xml"))
        topo = list(chk.parent.glob("*_analysis_topology.pdb")) or list((chk.parent / "assets").glob("*_start.pdb"))
        if not assets or not topo:
            return None
        system = mm.XmlSerializer.deserialize(assets[0].read_text())
        pdb = app.PDBFile(str(topo[0]))
        integ = mm.LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, DEFAULT_TIMESTEP_FS * unit.femtoseconds)
        ctx = mm.Context(system, integ, mm.Platform.getPlatformByName("Reference"))
        with open(chk, "rb") as fh:
            ctx.loadCheckpoint(fh.read())
        return int(ctx.getStepCount())
    except Exception:
        return None


def audit_run(rundir: Path, phase: str, *, timestep_fs: float, check_ckpt: bool) -> RunAudit:
    run_jsons = []
    for jf in sorted(rundir.glob("*.json")):
        try:
            d = json.loads(jf.read_text())
        except Exception:
            continue
        if "md_production_steps" in d:
            run_jsons.append((jf, d))

    primary = run_jsons[0][1] if run_jsons else {}
    json_steps = primary.get("md_production_steps_completed")
    stored = primary.get("checkpoint_path", "") or ""
    stored_ok = (REPO_ROOT / stored).exists() if stored else None

    state_steps = _last_state_step(rundir)
    chk = next(iter(rundir.glob("*_md.chk")), None)
    ckpt_steps = _checkpoint_step(chk) if (check_ckpt and chk) else None

    ns_done = round(state_steps * timestep_fs / 1e6, 1) if state_steps is not None else None
    return RunAudit(
        run=rel(rundir),
        phase=phase,
        json_files=len(run_jsons),
        json_steps=json_steps,
        state_steps=state_steps,
        ckpt_steps=ckpt_steps,
        ns_done=ns_done,
        json_matches_state=(json_steps == state_steps) if (json_steps is not None and state_steps is not None) else None,
        ckpt_matches_state=(ckpt_steps == state_steps) if (ckpt_steps is not None and state_steps is not None) else None,
        stored_path_ok=stored_ok,
        has_chk=chk is not None,
        has_system_xml=bool(list((rundir / "assets").glob("*_system.xml"))),
        status=primary.get("status"),
        diagnostic=(rundir / "DIAGNOSTIC.md").exists() or (rundir.parent / "DIAGNOSTIC.md").exists(),
    )


def collect(check_ckpt: bool, timestep_fs: float) -> list[RunAudit]:
    holo = [p for p in MD_RUNS.glob("*/rep_*") if p.is_dir() and p.parent.name != "apo"]
    apo = [p for p in APO_MD_RUNS.glob("*/rep_*") if p.is_dir()]
    out = []
    for d in sorted(holo):
        out.append(audit_run(d, "holo", timestep_fs=timestep_fs, check_ckpt=check_ckpt))
    for d in sorted(apo):
        out.append(audit_run(d, "apo", timestep_fs=timestep_fs, check_ckpt=check_ckpt))
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Read-only MD metadata consistency audit.")
    ap.add_argument("--check-checkpoints", action="store_true",
                    help="Also read each checkpoint's currentStep (needs OpenMM; slow).")
    ap.add_argument("--timestep-fs", type=float, default=DEFAULT_TIMESTEP_FS)
    ap.add_argument("--format", choices=["table", "csv"], default="table")
    ap.add_argument("--only-issues", action="store_true", help="Print only runs with inconsistencies.")
    args = ap.parse_args(argv)

    rows = collect(args.check_checkpoints, args.timestep_fs)
    shown = [r for r in rows if (not args.only_issues or not r.clean)]

    if args.format == "csv":
        w = csv.DictWriter(sys.stdout, fieldnames=list(asdict(rows[0]).keys()) if rows else [])
        w.writeheader()
        for r in shown:
            w.writerow(asdict(r))
        return 0

    hdr = f"{'run':<30} {'phase':<5} {'#json':>5} {'json':>10} {'state':>10} {'ckpt':>10} {'ns':>6} {'j=s':>4} {'c=s':>4} {'path':>5} {'chk':>4}"
    print(hdr)
    print("-" * len(hdr))
    for r in shown:
        print(f"{r.run.replace('results/md_runs/',''):<30} {r.phase:<5} {r.json_files:>5} "
              f"{str(r.json_steps):>10} {str(r.state_steps):>10} {str(r.ckpt_steps):>10} "
              f"{str(r.ns_done):>6} {str(r.json_matches_state):>4} {str(r.ckpt_matches_state):>4} "
              f"{str(r.stored_path_ok):>5} {str(r.has_chk):>4}")

    bad_path = sum(1 for r in rows if r.stored_path_ok is False)
    bad_json = sum(1 for r in rows if r.json_matches_state is False)
    bad_ckpt = sum(1 for r in rows if r.ckpt_matches_state is False)
    dup = sum(1 for r in rows if r.json_files != 1)
    print("-" * len(hdr))
    print(f"runs: {len(rows)} | stale paths: {bad_path} | json!=state: {bad_json} "
          f"| ckpt!=state: {bad_ckpt} | !=1 json/dir: {dup} | clean: {sum(r.clean for r in rows)}")
    for ph in ("holo", "apo"):
        ns = sorted({r.ns_done for r in rows if r.phase == ph and r.ns_done is not None})
        print(f"  {ph}: distinct ns on disk = {ns}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
