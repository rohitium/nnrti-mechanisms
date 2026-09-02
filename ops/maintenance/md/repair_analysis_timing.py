#!/usr/bin/env python3
"""Archive stale MD JSON/DCDs and write a 100 ns-corrected version.

The analysis DCD is the ground truth (see src/nnrti/analysis/md_timing.py). JSON and
state.csv often froze at 50–70 ns because the worker wrote status=ok after each
SLURM slice. DCD headers are also broken (nsavc=1, DELTA=1.0 ps).

This script:
  1. Copies the live JSON + analysis DCD to results/md_runs/_archive/
  2. Patches the live DCD header: NSTEP=50e6, DELTA in CHARMM AKMA
     (dt_ps = 100000/(n_frames-1); DELTA = dt_ps / 0.04888821)
  3. Writes a new JSON with md_production_steps_completed=50e6
  4. Patches aligned 4NCG DCDs the same way (header only; not re-archived)

State CSVs are left alone — they are incomplete energy logs, not fabricated.

    python3 ops/maintenance/md/repair_analysis_timing.py            # dry-run
    python3 ops/maintenance/md/repair_analysis_timing.py --apply
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import shutil
import struct
import sys
from pathlib import Path

REPO_ROOT = _repo_root()
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from nnrti.analysis.md_timing import DT_FS, infer_production_ns
from nnrti.paths import MANIFESTS, MD_RUNS, rel


def _repo_root() -> Path:
    """Nearest ancestor containing pyproject.toml."""
    for d in Path(__file__).resolve().parents:
        if (d / "pyproject.toml").is_file():
            return d
    raise RuntimeError("repository root not found from %s" % __file__)

ARCHIVE_DIR = MD_RUNS / "_archive"
ARCHIVE_LOG = MANIFESTS / "md_archive_log.csv"
REASON = "stale timing metadata + broken DCD header; replaced with 100 ns DCD-fingerprint version"
# CHARMM DCD DELTA is in AKMA time units. Readers (MDAnalysis) do dt_ps = DELTA * 0.04888821 * NSAVC.
AKMA_TO_PS = 0.04888821


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_header(path: Path) -> dict:
    raw = path.read_bytes()[:92]
    if raw[4:8] != b"CORD":
        raise ValueError(f"not a CORD DCD: {path}")
    nset = struct.unpack_from("<i", raw, 8)[0]
    nstart = struct.unpack_from("<i", raw, 12)[0]
    nsavc = struct.unpack_from("<i", raw, 16)[0]
    nstep = struct.unpack_from("<i", raw, 20)[0]
    delta = struct.unpack_from("<f", raw, 44)[0]
    return {"nset": nset, "nstart": nstart, "nsavc": nsavc, "nstep": nstep, "delta": delta}


def _patch_header(path: Path, *, nstep: int, delta_akma: float) -> None:
    with path.open("r+b") as handle:
        handle.seek(16)
        handle.write(struct.pack("<i", 1))  # NSAVC: DELTA is per-frame CHARMM/AKMA
        handle.seek(20)
        handle.write(struct.pack("<i", int(nstep)))
        handle.seek(44)
        handle.write(struct.pack("<f", float(delta_akma)))


def _append_archive_log(rows: list[dict]) -> None:
    exists = ARCHIVE_LOG.exists()
    ARCHIVE_LOG.parent.mkdir(parents=True, exist_ok=True)
    import csv

    with ARCHIVE_LOG.open("a", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["date", "reason", "original_path", "archived_path", "sha256", "bytes"],
        )
        if not exists:
            writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _holo_rep_dirs() -> list[Path]:
    return sorted(
        path
        for path in MD_RUNS.glob("*/rep_*")
        if path.is_dir() and path.parent.name not in {"apo", "_archive"} and "sherlock_rerun" not in path.parent.name
    )


def _primary_json(rundir: Path) -> Path | None:
    hits = []
    for path in sorted(rundir.glob("*.json")):
        if path.name.endswith("_quick.json"):
            continue
        try:
            payload = json.loads(path.read_text())
        except Exception:
            continue
        if "md_production_steps" in payload:
            hits.append(path)
    return hits[0] if hits else None


def _plan_one(rundir: Path) -> dict | None:
    json_path = _primary_json(rundir)
    if json_path is None:
        return None
    payload = json.loads(json_path.read_text())
    safe = str(payload.get("safe_label") or rundir.parent.name)
    rep = int(payload.get("replicate") or int(rundir.name.split("_")[-1]))
    dcd = rundir / f"{safe}_rep{rep:02d}_analysis.dcd"
    state_csv = rundir / f"{safe}_rep{rep:02d}_md_state.csv"
    aligned = rundir / f"{safe}_rep{rep:02d}_analysis_aligned_4ncg_ca.dcd"
    call = infer_production_ns(
        dcd_path=dcd if dcd.exists() else None,
        json_path=json_path,
        state_csv_path=state_csv if state_csv.exists() else None,
        mutation=str(payload.get("mutation") or safe),
        replicate=rep,
    )
    if call.n_frames is None or call.n_frames < 2:
        return None
    header = _read_header(dcd)
    production_ns = float(call.production_ns)
    production_steps = int(round(production_ns * 1e6 / DT_FS))
    delta_ps = (production_ns * 1000.0) / float(call.n_frames - 1)
    delta_akma = delta_ps / AKMA_TO_PS
    already = (
        int(payload.get("md_production_steps_completed") or 0) == production_steps
        and abs(float(header["delta"]) - delta_akma) < 1.0
        and int(header["nstep"]) == production_steps
        and bool(payload.get("timing_repaired"))
        and abs(float(payload.get("production_ns_completed") or 0) - production_ns) < 0.05
    )
    return {
        "rundir": rundir,
        "json_path": json_path,
        "payload": payload,
        "dcd": dcd,
        "aligned": aligned if aligned.exists() else None,
        "state_csv": state_csv if state_csv.exists() else None,
        "call": call,
        "old_header": header,
        "delta_ps": delta_ps,
        "delta_akma": delta_akma,
        "production_ns": production_ns,
        "production_steps": production_steps,
        "already": already,
    }


def _archive_copy(src: Path, apply: bool) -> tuple[Path, dict]:
    dest = ARCHIVE_DIR / src.relative_to(MD_RUNS)
    if apply:
        dest.parent.mkdir(parents=True, exist_ok=True)
        if not dest.exists():
            shutil.copy2(src, dest)
    return dest, {
        "date": _dt.date.today().isoformat(),
        "reason": REASON,
        "original_path": rel(src),
        "archived_path": rel(dest),
        "sha256": _sha256(dest if dest.exists() else src),
        "bytes": (dest if dest.exists() else src).stat().st_size,
    }


def _write_repaired_json(item: dict) -> None:
    payload = dict(item["payload"])
    heating = int(payload.get("md_heating_steps") or 0)
    call = item["call"]
    payload["md_production_steps"] = int(item["production_steps"])
    payload["md_production_steps_completed"] = int(item["production_steps"])
    payload["md_total_steps"] = heating + int(item["production_steps"])
    payload["production_ns_completed"] = float(item["production_ns"])
    payload["timing_repaired"] = _dt.date.today().isoformat()
    payload["timing_source"] = call.source
    payload["analysis_n_frames"] = call.n_frames
    payload["analysis_dt_ps"] = round(float(item["delta_ps"]), 6)
    payload["legacy_json_ns"] = call.json_ns
    payload["legacy_state_csv_ns"] = call.state_csv_ns
    payload["state_csv_stale"] = bool(
        call.state_csv_ns is not None and abs(float(call.state_csv_ns) - float(item["production_ns"])) > 5.0
    )
    item["json_path"].write_text(json.dumps(payload, indent=2) + "\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Archive stale MD timing files and write 100 ns-corrected JSON/DCD headers.")
    parser.add_argument("--apply", action="store_true", help="Perform changes (default: dry-run).")
    args = parser.parse_args(argv)

    items = [item for item in (_plan_one(path) for path in _holo_rep_dirs()) if item is not None]
    todo = [item for item in items if not item["already"]]
    mode = "APPLY" if args.apply else "DRY-RUN"
    print(f"=== analysis timing repair — {mode} ===")
    print(f"scanned {len(items)} holo reps; {len(todo)} need repair; {len(items) - len(todo)} already repaired\n")

    log_rows: list[dict] = []
    for item in todo:
        call = item["call"]
        old = item["old_header"]
        print(
            f"{rel(item['rundir']):<42} "
            f"{call.n_frames:>4}fr  json={call.json_ns or 0:5.1f}  csv={call.state_csv_ns or 0:5.1f}  "
            f"DELTA {old['delta']:.1f} -> {item['delta_akma']:.1f} AKMA ({item['delta_ps']:.2f} ps/frame)  "
            f"NSTEP {old['nstep']} -> {item['production_steps']}  ({item['production_ns']:.1f} ns)"
        )
        if not args.apply:
            continue
        dest_json = ARCHIVE_DIR / item["json_path"].relative_to(MD_RUNS)
        dest_dcd = ARCHIVE_DIR / item["dcd"].relative_to(MD_RUNS)
        if not dest_json.exists() or not dest_dcd.exists():
            _, json_log = _archive_copy(item["json_path"], apply=True)
            _, dcd_log = _archive_copy(item["dcd"], apply=True)
            log_rows.extend([json_log, dcd_log])
        _patch_header(item["dcd"], nstep=item["production_steps"], delta_akma=item["delta_akma"])
        if item["aligned"] is not None:
            _patch_header(item["aligned"], nstep=item["production_steps"], delta_akma=item["delta_akma"])
        _write_repaired_json(item)

    if args.apply and log_rows:
        _append_archive_log(log_rows)
        print(f"\narchived {len(log_rows)} files -> {rel(ARCHIVE_DIR)}")
        print(f"logged {len(log_rows)} rows -> {rel(ARCHIVE_LOG)}")

    if not args.apply and todo:
        print("\nRe-run with --apply to archive originals and write the corrected versions.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
