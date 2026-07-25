#!/usr/bin/env python3
"""Inventory local MD/FEP assets for the pmx NEQ panel."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.fep_jorgensen.mutations import MANUSCRIPT_PLANS, MutationLeg, unique_manuscript_legs
from scripts.fep_pmx.config import P0_LEGS


@dataclass(frozen=True)
class ReplicateAssets:
    leg_id: str
    phase: str
    genotype: str
    replicate: int
    start_pdb: str
    system_xml: str
    analysis_dcd: str
    md_final_pdb: str
    status: str
    trajectory_ready: bool


def _exists(path: Path) -> bool:
    return path.is_file()


def _holo_genotype(label: str) -> str:
    from scripts.fep_jorgensen.mutations import safe_label

    return safe_label(label)


def _replicate_row(leg: MutationLeg, replicate: int, phase: str) -> ReplicateAssets:
    mutation = leg.mutation
    if phase == "holo":
        genotype = _holo_genotype(leg.start_label)
        start = leg.input_complex_pdb(replicate)
        rep_dir = start.parent.parent
        prefix = "wt" if leg.start_label == "WT" else genotype
        dcd = rep_dir / f"{prefix}_rep{replicate:02d}_analysis.dcd"
        final_pdb = rep_dir / f"{prefix}_rep{replicate:02d}_md_final.pdb"
    else:
        from scripts.fep_jorgensen.mutations import apo_safe_label

        genotype = apo_safe_label(leg.start_label)
        start = leg.input_apo_pdb(replicate)
        rep_dir = start.parent.parent
        prefix = genotype
        dcd = rep_dir / f"{prefix}_rep{replicate:02d}_analysis.dcd"
        final_pdb = rep_dir / f"{prefix}_rep{replicate:02d}_md_final.pdb"

    json_candidates = sorted(rep_dir.glob(f"*{replicate:02d}*.json"))
    status = ""
    if json_candidates:
        try:
            payload = json.loads(json_candidates[0].read_text())
            status = str(payload.get("status", ""))
        except Exception:
            status = "json_error"

    system_xml = start.with_name(start.name.replace("_start.pdb", "_system.xml"))
    return ReplicateAssets(
        leg_id=leg.leg_id,
        phase=phase,
        genotype=genotype,
        replicate=replicate,
        start_pdb=str(start),
        system_xml=str(system_xml),
        analysis_dcd=str(dcd),
        md_final_pdb=str(final_pdb),
        status=status,
        trajectory_ready=_exists(start) and _exists(dcd),
    )


def collect_manifest(root: Path, *, replicates: tuple[int, ...] = (1, 2, 3)) -> list[ReplicateAssets]:
    rows: list[ReplicateAssets] = []
    for leg in unique_manuscript_legs():
        for phase in ("holo", "apo"):
            for rep in replicates:
                rows.append(_replicate_row(leg, rep, phase))
    return rows


def summarize(rows: list[ReplicateAssets]) -> dict[str, int]:
    return {
        "total_rows": len(rows),
        "trajectory_ready": sum(1 for row in rows if row.trajectory_ready),
        "start_pdb_only": sum(1 for row in rows if _exists(Path(row.start_pdb)) and not row.trajectory_ready),
        "missing_start": sum(1 for row in rows if not _exists(Path(row.start_pdb))),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Write pmx FEP asset manifest CSV/JSON.")
    parser.add_argument(
        "--root",
        type=Path,
        default=REPO_ROOT,
        help="Repository root (default: repo containing this script)",
    )
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=Path("results/analysis/fep_pmx/asset_manifest.csv"),
    )
    parser.add_argument(
        "--out-json",
        type=Path,
        default=Path("results/analysis/fep_pmx/asset_manifest_summary.json"),
    )
    args = parser.parse_args(argv)

    root = args.root.resolve()
    rows = collect_manifest(root)
    summary = summarize(rows)
    summary["p0_legs"] = list(P0_LEGS)
    summary["manuscript_targets"] = list(MANUSCRIPT_PLANS)

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))

    args.out_json.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Wrote {args.out_csv} ({len(rows)} rows)")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
