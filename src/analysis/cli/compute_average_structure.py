#!/usr/bin/env python3
"""Build a medoid structure from analysis trajectories for one mutation.

Example:
    conda run -n nnrti-prep python src/analysis/cli/compute_average_structure.py \
        --mutation WT \
        --window-ns 100 \
        --balance-replicates min_frames \
        --output-dir results/average_structures
"""
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass
class ReplicateMeta:
    mutation: str
    replicate: int
    output_json: Path
    topology_pdb: Path
    analysis_dcd: Path
    total_ns: float
    timing_source: str


def _steps_to_ns(steps: float | int | None, timestep_fs: float = 2.0) -> float:
    try:
        v = float(steps)
    except Exception:
        return np.nan
    if not np.isfinite(v):
        return np.nan
    return float(v * timestep_fs / 1_000_000.0)


def _resolve_local_path(path_like: str | Path | None, repo_root: Path) -> Path | None:
    if path_like is None:
        return None
    p = Path(str(path_like))
    if p.exists():
        return p
    marker = "nnrti-mechanisms/"
    text = str(p)
    if marker in text:
        mapped = repo_root / text.split(marker, 1)[1]
        if mapped.exists():
            return mapped
    rel = repo_root / str(p)
    if rel.exists():
        return rel
    return p


def _infer_total_ns_from_state_csv(path: Path | None) -> float:
    if path is None or not path.exists():
        return np.nan
    try:
        sdf = pd.read_csv(path)
    except Exception:
        return np.nan
    if sdf.empty:
        return np.nan
    step_col = None
    for c in ['#"Step"', '"#Step"', "Step"]:
        if c in sdf.columns:
            step_col = c
            break
    if step_col is None:
        return np.nan
    steps = pd.to_numeric(sdf[step_col], errors="coerce").dropna()
    if steps.empty:
        return np.nan
    return _steps_to_ns(float(steps.max()), timestep_fs=2.0)


def _load_replicate_meta(manifest_csv: Path, mutation: str) -> list[ReplicateMeta]:
    repo_root = Path(__file__).resolve().parents[3]
    mf = pd.read_csv(manifest_csv)
    req_cols = {"mutation", "replicate", "output_json"}
    missing = req_cols - set(mf.columns)
    if missing:
        raise ValueError(f"Manifest missing required columns: {sorted(missing)}")

    mutation_key = str(mutation).strip().upper()
    out: list[ReplicateMeta] = []
    rows = mf[mf["mutation"].astype(str).str.upper() == mutation_key].copy()

    for _, row in rows.sort_values(["mutation", "replicate"]).iterrows():
        rep = int(pd.to_numeric(row["replicate"], errors="coerce"))
        out_json = _resolve_local_path(row["output_json"], repo_root=repo_root)
        if out_json is None or not out_json.exists():
            continue
        try:
            data = json.loads(out_json.read_text())
        except Exception:
            continue

        topo = _resolve_local_path(data.get("analysis_topology_pdb"), repo_root=repo_root)
        dcd = _resolve_local_path(data.get("analysis_dcd"), repo_root=repo_root)
        if topo is None or dcd is None or (not topo.exists()) or (not dcd.exists()):
            continue

        ns_json = _steps_to_ns(
            data.get("md_production_steps_completed", data.get("md_production_steps")),
            timestep_fs=2.0,
        )
        state_csv = _resolve_local_path(data.get("state_csv"), repo_root=repo_root)
        ns_state = _infer_total_ns_from_state_csv(state_csv)
        has_state = bool(np.isfinite(ns_state) and ns_state > 0)
        has_json = bool(np.isfinite(ns_json) and ns_json > 0)
        if has_state and has_json:
            total_ns = float(max(ns_state, ns_json))
            timing_source = "state_csv" if ns_state >= ns_json else "json_steps_gt_state_csv"
        elif has_state:
            total_ns = float(ns_state)
            timing_source = "state_csv"
        elif has_json:
            total_ns = float(ns_json)
            timing_source = "json_steps"
        else:
            total_ns = np.nan
            timing_source = "unknown"

        out.append(
            ReplicateMeta(
                mutation=str(row["mutation"]),
                replicate=rep,
                output_json=out_json,
                topology_pdb=topo,
                analysis_dcd=dcd,
                total_ns=total_ns,
                timing_source=timing_source,
            )
        )
    return sorted(out, key=lambda m: m.replicate)


def _atom_signature(topology) -> tuple[tuple[int, int, str, str, str], ...]:
    sig: list[tuple[int, int, str, str, str]] = []
    for atom in topology.atoms:
        res = atom.residue
        elem = str(getattr(atom.element, "symbol", "") or "")
        sig.append((int(res.chain.index), int(res.resSeq), str(res.name), str(atom.name), elem))
    return tuple(sig)


def _make_molecules_whole_nojump(traj) -> None:
    traj.make_molecules_whole(inplace=True)
    if getattr(traj, "unitcell_lengths", None) is None:
        return
    mol_indices = [
        np.asarray([atom.index for atom in mol], dtype=int)
        for mol in traj.topology.find_molecules()
    ]
    for frame_i in range(1, traj.n_frames):
        box = traj.unitcell_lengths[frame_i]
        if box is None or not np.all(np.isfinite(box)) or not np.all(box > 0):
            continue
        for mol_idx in mol_indices:
            if mol_idx.size == 0:
                continue
            prev_com = traj.xyz[frame_i - 1, mol_idx].mean(axis=0)
            curr_com = traj.xyz[frame_i, mol_idx].mean(axis=0)
            shift = -box * np.round((curr_com - prev_com) / box)
            traj.xyz[frame_i, mol_idx] += shift


def _window_trajectory(traj, total_ns: float, window_ns: float):
    n_frames = int(traj.n_frames)
    t_ns = np.linspace(0.0, float(total_ns), n_frames)
    keep = t_ns <= float(window_ns)
    if int(np.sum(keep)) < 2:
        keep = np.ones(n_frames, dtype=bool)
    keep_idx = np.where(keep)[0].astype(int)
    return traj[keep], t_ns[keep], keep_idx, int(n_frames), int(np.sum(keep))


def _uniform_subsample_indices(n_frames: int, n_keep: int) -> np.ndarray:
    n_frames = int(n_frames)
    n_keep = int(n_keep)
    if n_keep <= 0:
        raise ValueError("n_keep must be > 0")
    if n_keep >= n_frames:
        return np.arange(n_frames, dtype=int)
    idx = np.floor(np.arange(n_keep, dtype=float) * float(n_frames) / float(n_keep)).astype(int)
    idx = np.clip(idx, 0, n_frames - 1)
    return idx


def _default_output_name(mutation: str) -> str:
    return f"{str(mutation).strip()}_medoid_structure.pdb"


def _mutation_dir_name(mutation: str) -> str:
    return str(mutation).strip().replace("+", "_")


def main() -> int:
    parser = argparse.ArgumentParser(description="Compute a medoid structure for one mutation.")
    parser.add_argument("--manifest", type=Path, default=Path("manifests/md_manifest.csv"))
    parser.add_argument("--mutation", type=str, default="WT")
    parser.add_argument("--window-ns", type=float, default=100.0)
    parser.add_argument(
        "--align-selection",
        type=str,
        default="protein and backbone and name CA",
        help="mdtraj DSL selection used for frame superposition.",
    )
    parser.add_argument(
        "--medoid-selection",
        type=str,
        default="",
        help="mdtraj DSL selection used for medoid distance metric (default: --align-selection).",
    )
    parser.add_argument(
        "--balance-replicates",
        choices=["min_frames", "none", "fixed"],
        default="min_frames",
        help=(
            "How to prevent long replicates from dominating medoid selection: "
            "'min_frames' downsamples every replicate to shortest frame count "
            "after windowing; 'fixed' uses --fixed-frames-per-rep; 'none' keeps all frames."
        ),
    )
    parser.add_argument(
        "--fixed-frames-per-rep",
        type=int,
        default=0,
        help="Required when --balance-replicates fixed.",
    )
    parser.add_argument("--output-dir", type=Path, default=Path("results/average_structures"))
    parser.add_argument(
        "--flat-output",
        action="store_true",
        help="Write files directly into --output-dir instead of mutation subfolder.",
    )
    parser.add_argument("--output-name", type=str, default="")
    args = parser.parse_args()

    import mdtraj as md

    metas = _load_replicate_meta(args.manifest, mutation=args.mutation)
    if not metas:
        raise ValueError(f"No usable replicates found for mutation {args.mutation!r}")

    out_dir = args.output_dir if args.flat_output else (args.output_dir / _mutation_dir_name(args.mutation))
    out_dir.mkdir(parents=True, exist_ok=True)
    out_name = args.output_name.strip() or _default_output_name(args.mutation)
    out_pdb = out_dir / out_name

    ref_frame = None
    ref_align_idx = None
    ref_sig = None
    total_frames = 0
    joined_trajs: list[object] = []
    frame_rows: list[dict[str, object]] = []
    prepped: list[dict[str, object]] = []
    rep_rows: list[dict[str, object]] = []

    for m in metas:
        traj = md.load_dcd(str(m.analysis_dcd), top=str(m.topology_pdb))
        usable_ns = float(m.total_ns) if np.isfinite(m.total_ns) and m.total_ns > 0 else float(args.window_ns)
        traj, t_sel, window_frame_idx, n_total, n_window = _window_trajectory(
            traj,
            total_ns=usable_ns,
            window_ns=float(args.window_ns),
        )
        _make_molecules_whole_nojump(traj)
        prepped.append(
            {
                "meta": m,
                "traj": traj,
                "n_total": int(n_total),
                "n_window": int(n_window),
                "usable_ns": float(usable_ns),
                "time_ns": np.asarray(t_sel, dtype=float),
                "window_frame_idx": np.asarray(window_frame_idx, dtype=int),
            }
        )

    if not prepped:
        raise ValueError("No trajectories available after preprocessing.")

    window_counts = [int(entry["traj"].n_frames) for entry in prepped]
    if args.balance_replicates == "min_frames":
        target_per_rep = int(min(window_counts))
    elif args.balance_replicates == "fixed":
        if int(args.fixed_frames_per_rep) <= 0:
            raise ValueError("--fixed-frames-per-rep must be > 0 when --balance-replicates fixed")
        target_per_rep = int(args.fixed_frames_per_rep)
        too_short = [
            f"{entry['meta'].mutation} rep{entry['meta'].replicate} ({entry['traj'].n_frames} frames)"
            for entry in prepped
            if entry["traj"].n_frames < target_per_rep
        ]
        if too_short:
            raise ValueError(
                "Requested fixed frame count exceeds available frames for: " + ", ".join(too_short)
            )
    else:
        target_per_rep = -1

    global_frame = 0
    for entry in prepped:
        m = entry["meta"]
        traj = entry["traj"]
        n_total = int(entry["n_total"])
        n_window = int(entry["n_window"])
        usable_ns = float(entry["usable_ns"])
        t_sel = np.asarray(entry["time_ns"], dtype=float)
        window_frame_idx = np.asarray(entry["window_frame_idx"], dtype=int)

        if args.balance_replicates in {"min_frames", "fixed"}:
            keep_idx = _uniform_subsample_indices(int(traj.n_frames), int(target_per_rep))
            traj = traj[keep_idx]
            t_sel = t_sel[keep_idx]
            window_frame_idx = window_frame_idx[keep_idx]
        else:
            keep_idx = np.arange(int(traj.n_frames), dtype=int)
        n_balanced = int(traj.n_frames)

        align_idx = traj.topology.select(str(args.align_selection))
        if align_idx.size < 3:
            raise ValueError(
                f"Alignment selection '{args.align_selection}' produced too few atoms "
                f"for {m.mutation} rep{m.replicate}."
            )

        sig = _atom_signature(traj.topology)
        if ref_sig is None:
            ref_sig = sig
            ref_frame = traj[0]
            ref_align_idx = align_idx
        else:
            if sig != ref_sig:
                raise ValueError(
                    f"Topology atom order mismatch for {m.mutation} rep{m.replicate}; "
                    "cannot compare non-identical atom orderings."
                )
            traj.superpose(ref_frame, atom_indices=align_idx, ref_atom_indices=ref_align_idx)

        joined_trajs.append(traj)
        for local_idx in range(traj.n_frames):
            frame_rows.append(
                {
                    "global_frame_index": int(global_frame),
                    "mutation": m.mutation,
                    "replicate": int(m.replicate),
                    "frame_local_index": int(local_idx),
                    "frame_window_index": int(window_frame_idx[local_idx]),
                    "time_ns": float(t_sel[local_idx]),
                }
            )
            global_frame += 1
        total_frames += int(traj.n_frames)

        rep_rows.append(
            {
                "mutation": m.mutation,
                "replicate": int(m.replicate),
                "analysis_dcd": str(m.analysis_dcd),
                "analysis_topology_pdb": str(m.topology_pdb),
                "n_frames_total": int(n_total),
                "n_frames_window": int(n_window),
                "n_frames_used": int(n_balanced),
                "total_ns_used": float(usable_ns),
                "timing_source": m.timing_source,
                "balance_mode": str(args.balance_replicates),
                "target_frames_per_rep": int(target_per_rep) if target_per_rep > 0 else "",
                "first_time_ns": float(t_sel[0]) if t_sel.size else np.nan,
                "last_time_ns": float(t_sel[-1]) if t_sel.size else np.nan,
                "downsample_keep_index_first": int(keep_idx[0]) if keep_idx.size else -1,
                "downsample_keep_index_last": int(keep_idx[-1]) if keep_idx.size else -1,
            }
        )

    if total_frames < 1 or ref_frame is None or not joined_trajs:
        raise ValueError("No frames available for medoid construction.")

    medoid_sel = str(args.medoid_selection).strip() or str(args.align_selection)
    full = md.join(joined_trajs, check_topology=True)
    medoid_idx = full.topology.select(medoid_sel)
    if medoid_idx.size < 3:
        raise ValueError(
            f"Medoid selection '{medoid_sel}' produced too few atoms ({medoid_idx.size})."
        )

    dist_sums = np.zeros(full.n_frames, dtype=np.float64)
    for i in range(full.n_frames):
        d = md.rmsd(full, full, i, atom_indices=medoid_idx)
        dist_sums[i] = float(np.sum(d))
    best_idx = int(np.argmin(dist_sums))
    med = full[best_idx]
    med.save_pdb(str(out_pdb))

    frame_df = pd.DataFrame(frame_rows).sort_values("global_frame_index").reset_index(drop=True)
    best = frame_df.iloc[best_idx].to_dict()

    rep_csv = out_dir / f"{out_pdb.stem}_replicate_audit.csv"
    pd.DataFrame(rep_rows).sort_values(["replicate"]).to_csv(rep_csv, index=False)
    frame_csv = out_dir / f"{out_pdb.stem}_frame_index.csv"
    frame_df.to_csv(frame_csv, index=False)

    summary_json = out_dir / f"{out_pdb.stem}_summary.json"
    summary_json.write_text(
        json.dumps(
            {
                "mutation": str(args.mutation),
                "window_ns": float(args.window_ns),
                "alignment_selection": str(args.align_selection),
                "medoid_selection": medoid_sel,
                "balance_replicates": str(args.balance_replicates),
                "target_frames_per_rep": int(target_per_rep) if target_per_rep > 0 else None,
                "n_replicates": int(len(metas)),
                "n_frames_used": int(total_frames),
                "medoid_global_frame_index": int(best_idx),
                "medoid_rmsd_sum_nm": float(dist_sums[best_idx]),
                "medoid_replicate": int(best["replicate"]),
                "medoid_frame_local_index": int(best["frame_local_index"]),
                "medoid_frame_window_index": int(best["frame_window_index"]),
                "medoid_time_ns": float(best["time_ns"]),
                "output_pdb": str(out_pdb),
                "replicate_audit_csv": str(rep_csv),
                "frame_index_csv": str(frame_csv),
            },
            indent=2,
        )
    )

    print(f"Saved {out_pdb}")
    print(f"Saved {rep_csv}")
    print(f"Saved {frame_csv}")
    print(f"Saved {summary_json}")
    print(f"Frames considered for medoid: {total_frames}")
    print(
        "Medoid frame: "
        f"rep{int(best['replicate'])} "
        f"local={int(best['frame_local_index'])} "
        f"time_ns={float(best['time_ns']):.3f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
