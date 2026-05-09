#!/usr/bin/env python3
"""Plot a 100 ns residue-to-DOR distance triplet story on a common time grid."""
from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..susceptibility import load_dor_susceptibilities


@dataclass
class ReplicateMeta:
    mutation: str
    replicate: int
    output_json: Path
    topology_pdb: Path
    analysis_dcd: Path
    total_ns: float
    timing_source: str


MUTATION_COLORS = {
    "WT": "#333333",
    "Y181C": "#4c78a8",
    "Y188L": "#e45756",
    "V106I": "#4c78a8",
    "V106A": "#e45756",
    "G190A": "#4c78a8",
    "G190E": "#e45756",
}


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot residue-to-DOR triplet stories over 100 ns.")
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("manifests/md_manifest.csv"),
    )
    parser.add_argument(
        "--susceptibility-xlsx",
        type=Path,
        default=Path("data/DRM-susceptibilities.csv.xlsx"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/analysis/triplet_geometry_story_100ns_v106"),
    )
    parser.add_argument("--triplet", type=str, default="WT,V106I,V106A")
    parser.add_argument("--auth-resseq", type=int, default=105)
    parser.add_argument("--metric-name", type=str, default="SER105-DOR")
    parser.add_argument("--ylabel", type=str, default="Min SER105-DOR Distance (Å)")
    parser.add_argument("--output-prefix", type=str, default="triplet_story_100ns_WT_V106I_V106A_SER105")
    parser.add_argument("--max-time-ns", type=float, default=100.0)
    parser.add_argument(
        "--force-total-ns",
        type=float,
        default=None,
        help="If set, ignore per-replicate timing metadata and map each trajectory uniformly onto this duration.",
    )
    parser.add_argument("--resid-offset", type=int, default=-3)
    parser.add_argument("--ligand-resname", type=str, default="2KW")
    parser.add_argument(
        "--triplet-colors",
        type=str,
        default="",
        help="Optional comma-separated colors to apply to the triplet in order.",
    )
    return parser.parse_args()


def _load_fold_map(xlsx_path: Path) -> dict[str, float]:
    df = load_dor_susceptibilities(xlsx_path)
    fold_map = {str(row["mutation"]): float(row["dor_fold_reduction"]) for _, row in df.iterrows()}
    fold_map["WT"] = 1.0
    return fold_map


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


def _resolve_dcd_with_fallback(dcd_path: Path | None) -> Path | None:
    if dcd_path is None:
        return None
    if dcd_path.exists():
        return dcd_path
    name = dcd_path.name
    candidates: list[Path] = []
    if name.endswith("_analysis.dcd"):
        candidates.append(dcd_path.with_name(name.replace("_analysis.dcd", "_analysis.10ns.bak")))
    candidates.append(dcd_path.with_suffix(dcd_path.suffix + ".bak"))
    candidates.append(dcd_path.with_name(name + ".bak"))
    for c in candidates:
        if c.exists():
            return c
    return dcd_path


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


def _load_replicate_meta(manifest_csv: Path, needed_mutations: set[str]) -> list[ReplicateMeta]:
    repo_root = Path(__file__).resolve().parents[3]
    mf = pd.read_csv(manifest_csv)
    req_cols = {"mutation", "replicate", "output_json"}
    missing = req_cols - set(mf.columns)
    if missing:
        raise ValueError(f"Manifest missing required columns: {sorted(missing)}")

    out: list[ReplicateMeta] = []
    for _, row in mf.sort_values(["mutation", "replicate"]).iterrows():
        mutation = str(row["mutation"])
        if mutation not in needed_mutations:
            continue
        replicate = int(pd.to_numeric(row["replicate"], errors="coerce"))
        out_json = _resolve_local_path(row["output_json"], repo_root=repo_root)
        if out_json is None or not out_json.exists():
            continue
        try:
            data = json.loads(out_json.read_text())
        except Exception:
            continue

        topo = _resolve_local_path(data.get("analysis_topology_pdb"), repo_root=repo_root)
        dcd = _resolve_local_path(data.get("analysis_dcd"), repo_root=repo_root)
        dcd = _resolve_dcd_with_fallback(dcd)
        if topo is None or dcd is None or (not topo.exists()) or (not dcd.exists()):
            continue

        ns_json = _steps_to_ns(data.get("md_production_steps_completed", data.get("md_production_steps")), timestep_fs=2.0)
        state_csv = _resolve_local_path(data.get("state_csv"), repo_root=repo_root)
        ns_state = _infer_total_ns_from_state_csv(state_csv)

        has_state = bool(np.isfinite(ns_state) and ns_state > 0)
        has_json = bool(np.isfinite(ns_json) and ns_json > 0)
        if has_state and has_json:
            if float(ns_state) >= float(ns_json):
                total_ns = float(ns_state)
                timing_source = "state_csv"
            else:
                total_ns = float(ns_json)
                timing_source = "json_steps_gt_state_csv"
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
                mutation=mutation,
                replicate=replicate,
                output_json=out_json,
                topology_pdb=topo,
                analysis_dcd=dcd,
                total_ns=total_ns,
                timing_source=timing_source,
            )
        )
    return out


def _largest_protein_chain(topology):
    protein_chains = [c for c in topology.chains if sum(r.is_protein for r in c.residues) > 0]
    if not protein_chains:
        raise ValueError("No protein chain found in topology")
    return max(protein_chains, key=lambda c: sum(r.is_protein for r in c.residues))


def _protein_residue_index(topology, chain_idx: int, resseq: int) -> int:
    matches = [
        int(residue.index)
        for residue in topology.residues
        if residue.is_protein and int(residue.chain.index) == int(chain_idx) and int(residue.resSeq) == int(resseq)
    ]
    if len(matches) != 1:
        raise ValueError(f"Could not map protein residue {resseq} on chain index {chain_idx}")
    return int(matches[0])


def _ligand_residue_index(topology, ligand_resname: str) -> int:
    matches = [int(residue.index) for residue in topology.residues if str(residue.name).strip() == str(ligand_resname).strip()]
    if len(matches) != 1:
        raise ValueError(f"Could not map ligand residue {ligand_resname}")
    return int(matches[0])


def _build_common_time_grid(trace_df: pd.DataFrame, max_time_ns: float) -> np.ndarray:
    times = pd.to_numeric(trace_df["time_ns"], errors="coerce").dropna().to_numpy(dtype=float)
    times = np.unique(np.round(times, 6))
    times = times[(times >= 0.0) & (times <= float(max_time_ns))]
    if times.size == 0:
        return np.linspace(0.0, float(max_time_ns), num=101, dtype=float)
    if times[0] > 0.0:
        times = np.concatenate([np.array([0.0], dtype=float), times])
    if times[-1] < float(max_time_ns):
        times = np.concatenate([times, np.array([float(max_time_ns)], dtype=float)])
    return times


def _compute_residue_trace(
    analysis_dcd: str | Path,
    topology_pdb: str | Path,
    total_ns: float,
    window_ns: float,
    *,
    auth_resseq: int,
    resid_offset: int,
    ligand_resname: str,
) -> tuple[np.ndarray, np.ndarray]:
    try:
        return _compute_residue_trace_mdtraj(
            analysis_dcd=analysis_dcd,
            topology_pdb=topology_pdb,
            total_ns=total_ns,
            window_ns=window_ns,
            auth_resseq=auth_resseq,
            resid_offset=resid_offset,
            ligand_resname=ligand_resname,
        )
    except ModuleNotFoundError as exc:
        if exc.name != "mdtraj":
            raise
        return _compute_residue_trace_mdanalysis(
            analysis_dcd=analysis_dcd,
            topology_pdb=topology_pdb,
            total_ns=total_ns,
            window_ns=window_ns,
            auth_resseq=auth_resseq,
            resid_offset=resid_offset,
            ligand_resname=ligand_resname,
        )


def _compute_residue_trace_mdtraj(
    analysis_dcd: str | Path,
    topology_pdb: str | Path,
    total_ns: float,
    window_ns: float,
    *,
    auth_resseq: int,
    resid_offset: int,
    ligand_resname: str,
) -> tuple[np.ndarray, np.ndarray]:
    import mdtraj as md

    traj = md.load_dcd(str(analysis_dcd), top=str(topology_pdb))
    n_frames = int(traj.n_frames)
    if n_frames < 2:
        raise ValueError(f"Too few frames in {analysis_dcd}")

    total_ns = float(total_ns) if np.isfinite(total_ns) and total_ns > 0 else float(window_ns)
    t_ns = np.linspace(0.0, total_ns, n_frames)
    keep = t_ns <= float(window_ns)
    if int(np.sum(keep)) < 2:
        keep = np.ones(n_frames, dtype=bool)
    sub = traj[keep]
    t_sel = t_ns[keep]
    chain = _largest_protein_chain(sub.topology)
    protein_res_idx = _protein_residue_index(sub.topology, int(chain.index), int(auth_resseq) + int(resid_offset))
    ligand_res_idx = _ligand_residue_index(sub.topology, ligand_resname)
    distances_nm, _ = md.compute_contacts(
        sub,
        contacts=np.asarray([[protein_res_idx, ligand_res_idx]], dtype=int),
        scheme="closest-heavy",
        periodic=True,
    )
    return t_sel.astype(float), (distances_nm[:, 0] * 10.0).astype(float)


def _compute_residue_trace_mdanalysis(
    analysis_dcd: str | Path,
    topology_pdb: str | Path,
    total_ns: float,
    window_ns: float,
    *,
    auth_resseq: int,
    resid_offset: int,
    ligand_resname: str,
) -> tuple[np.ndarray, np.ndarray]:
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import distance_array

    u = mda.Universe(str(topology_pdb), str(analysis_dcd))
    n_frames = len(u.trajectory)
    if n_frames < 2:
        raise ValueError(f"Too few frames in {analysis_dcd}")

    total_ns = float(total_ns) if np.isfinite(total_ns) and total_ns > 0 else float(window_ns)
    t_ns = np.linspace(0.0, total_ns, n_frames)
    target_resid = int(auth_resseq) + int(resid_offset)

    protein = u.select_atoms("protein")
    segments = sorted(protein.segments, key=lambda seg: len(seg.atoms), reverse=True)
    if not segments:
        raise ValueError("No protein segment found in topology")
    segid = str(segments[0].segid)
    residue = u.select_atoms(f"protein and segid {segid} and resid {target_resid} and not name H*")
    ligand = u.select_atoms(f"resname {ligand_resname} and not name H*")
    if len(residue) == 0:
        raise ValueError(f"Could not map protein residue {target_resid} on segment {segid}")
    if len(ligand) == 0:
        raise ValueError(f"Could not map ligand residue {ligand_resname}")

    times: list[float] = []
    distances: list[float] = []
    for frame_idx, ts in enumerate(u.trajectory):
        time_ns = float(t_ns[frame_idx])
        if time_ns > float(window_ns):
            continue
        d = distance_array(residue.positions, ligand.positions, box=ts.dimensions)
        times.append(time_ns)
        distances.append(float(np.nanmin(d)))
    if len(times) < 2:
        raise ValueError(f"Too few frames retained in {analysis_dcd}")
    return np.asarray(times, dtype=float), np.asarray(distances, dtype=float)


def _interpolate_replicates_to_grid(trace_df: pd.DataFrame, metric: str, time_grid: np.ndarray) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for (mutation, replicate), rep_df in trace_df.groupby(["mutation", "replicate"], sort=True):
        rep_df = rep_df.sort_values("time_ns", kind="stable").copy()
        x = pd.to_numeric(rep_df["time_ns"], errors="coerce").to_numpy(dtype=float)
        y = pd.to_numeric(rep_df[metric], errors="coerce").to_numpy(dtype=float)
        ok = np.isfinite(x) & np.isfinite(y)
        x = x[ok]
        y = y[ok]
        if x.size == 0:
            continue
        x_unique, unique_idx = np.unique(x, return_index=True)
        y_unique = y[unique_idx]
        if x_unique.size < 2:
            continue
        y_interp = np.interp(time_grid, x_unique, y_unique, left=np.nan, right=np.nan)
        y_interp[(time_grid < x_unique.min()) | (time_grid > x_unique.max())] = np.nan
        rows.append(
            pd.DataFrame(
                {
                    "mutation": str(mutation),
                    "replicate": int(replicate),
                    "time_ns": time_grid.astype(float),
                    metric: y_interp.astype(float),
                }
            )
        )
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["mutation", "replicate", "time_ns", metric])


def main() -> int:
    args = _parse_args()
    if not args.manifest.exists():
        raise FileNotFoundError(args.manifest)
    if not args.susceptibility_xlsx.exists():
        raise FileNotFoundError(args.susceptibility_xlsx)

    triplet = [token.strip() for token in str(args.triplet).split(",") if token.strip()]
    if len(triplet) != 3:
        raise ValueError("--triplet must contain exactly three comma-separated mutations")
    triplet_colors = [token.strip() for token in str(args.triplet_colors).split(",") if token.strip()]
    if triplet_colors and len(triplet_colors) != len(triplet):
        raise ValueError("--triplet-colors must provide exactly one color per triplet entry")

    out_tables = args.output_dir / "tables"
    out_plots = args.output_dir / "plots"
    out_config = args.output_dir / "config"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_plots.mkdir(parents=True, exist_ok=True)
    out_config.mkdir(parents=True, exist_ok=True)

    fold_map = _load_fold_map(args.susceptibility_xlsx)
    metas = _load_replicate_meta(args.manifest, needed_mutations=set(triplet))
    if not metas:
        raise ValueError("No valid replicate metadata found for requested triplet")

    raw_metric_frames: list[pd.DataFrame] = []
    interp_metric_frames: list[pd.DataFrame] = []
    mean_rows: list[pd.DataFrame] = []

    metric_rows: list[pd.DataFrame] = []
    for meta in sorted(metas, key=lambda x: (x.mutation, x.replicate)):
        t_sel, trace = _compute_residue_trace(
            analysis_dcd=str(meta.analysis_dcd),
            topology_pdb=str(meta.topology_pdb),
            total_ns=float(args.force_total_ns) if args.force_total_ns is not None else float(meta.total_ns),
            window_ns=float(args.max_time_ns),
            auth_resseq=int(args.auth_resseq),
            resid_offset=int(args.resid_offset),
            ligand_resname=str(args.ligand_resname),
        )
        sub = pd.DataFrame(
            {
                "mutation": str(meta.mutation),
                "replicate": int(meta.replicate),
                "time_ns": t_sel.astype(float),
            }
        )
        sub["metric"] = str(args.metric_name)
        sub["metric_value"] = trace
        metric_rows.append(sub)
    metric_df = pd.concat(metric_rows, ignore_index=True)
    metric_df = metric_df.sort_values(["mutation", "replicate", "time_ns"], kind="stable").reset_index(drop=True)
    metric_df.to_csv(out_tables / "trace_values_raw.csv", index=False)
    raw_metric_frames.append(metric_df)

    time_grid = _build_common_time_grid(metric_df, max_time_ns=float(args.max_time_ns))

    interp_df = _interpolate_replicates_to_grid(
        metric_df.rename(columns={"metric_value": str(args.metric_name)})[
            ["mutation", "replicate", "time_ns", str(args.metric_name)]
        ],
        str(args.metric_name),
        time_grid,
    )
    interp_df["metric"] = str(args.metric_name)
    interp_metric_frames.append(interp_df.rename(columns={str(args.metric_name): "metric_value"}))

    mean_df = (
        interp_df.groupby(["mutation", "time_ns"], as_index=False)[str(args.metric_name)]
        .agg(["mean", "std", "count"])
        .reset_index()
        .rename(columns={"mean": "metric_mean", "std": "metric_std", "count": "n_replicates"})
    )
    mean_df["metric"] = str(args.metric_name)
    mean_df["metric_sem"] = (
        mean_df["metric_std"].fillna(0.0) / mean_df["n_replicates"].clip(lower=1).pow(0.5)
    ).astype(float)
    mean_rows.append(mean_df[["mutation", "time_ns", "metric", "metric_mean", "metric_std", "metric_sem", "n_replicates"]])

    pd.concat(raw_metric_frames, ignore_index=True).to_csv(out_tables / "trace_values_raw_metrics.csv", index=False)
    pd.concat(interp_metric_frames, ignore_index=True).to_csv(out_tables / "trace_values.csv", index=False)
    mean_df_all = pd.concat(mean_rows, ignore_index=True)
    mean_df_all.to_csv(out_tables / "mean_traces.csv", index=False)

    fig, ax = plt.subplots(figsize=(15.8, 6.4), constrained_layout=True)

    xmax = 0.0
    metric_name = str(args.metric_name)
    ylabel = str(args.ylabel)
    sub = mean_df_all[mean_df_all["metric"].astype(str) == metric_name].copy()
    for idx, mutation in enumerate(triplet):
        color = triplet_colors[idx] if triplet_colors else MUTATION_COLORS.get(mutation, "#555555")
        mut_mean = sub[sub["mutation"].astype(str) == mutation].copy()
        x = mut_mean["time_ns"].to_numpy(dtype=float)
        y = mut_mean["metric_mean"].to_numpy(dtype=float)
        sem = mut_mean["metric_sem"].to_numpy(dtype=float)
        xmax = max(xmax, float(np.nanmax(x)) if len(x) else 0.0)
        fold = fold_map.get(mutation, float("nan"))
        label = f"{mutation} ({fold:.1f}x)" if pd.notna(fold) else str(mutation)
        ax.plot(x, y, color=color, linewidth=2.8, alpha=0.95, label=label)
        lo = y - sem
        hi = y + sem
        ok = np.isfinite(lo) & np.isfinite(hi)
        if np.any(ok):
            ax.fill_between(x[ok], lo[ok], hi[ok], color=color, alpha=0.16, linewidth=0)
    ax.axhline(4.0, color="#666666", linestyle=":", linewidth=1.8, label="Contact cutoff (4 Å)")
    ax.set_xlim(0.0, xmax if xmax > 0 else float(args.max_time_ns))
    ax.set_xlabel("Time (ns)", fontsize=18)
    ax.set_ylabel(ylabel, fontsize=18)
    ax.set_title(metric_name, fontsize=20, fontweight="bold")
    ax.tick_params(axis="both", labelsize=19)
    ax.grid(alpha=0.22, linestyle=":")
    ax.legend(loc="upper right", frameon=True, fontsize=14)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    png = out_plots / f"{str(args.output_prefix)}.png"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    (out_config / "run_config.json").write_text(
        json.dumps(
            {
                "manifest": str(args.manifest),
                "susceptibility_xlsx": str(args.susceptibility_xlsx),
                "output_dir": str(args.output_dir),
                "triplet": triplet,
                "max_time_ns": float(args.max_time_ns),
                "resid_offset": int(args.resid_offset),
                "ligand_resname": str(args.ligand_resname),
                "auth_resseq": int(args.auth_resseq),
                "metrics": [str(args.metric_name)],
                "ylabel": str(args.ylabel),
                "output_prefix": str(args.output_prefix),
                "triplet_colors": triplet_colors,
                "force_total_ns": None if args.force_total_ns is None else float(args.force_total_ns),
                "aggregation": "replicates interpolated onto common time grid before mean/SEM",
                "timing_mode": (
                    "forced_uniform_duration"
                    if args.force_total_ns is not None
                    else "per_replicate_metadata_duration"
                ),
                "common_time_grid_size": int(len(time_grid)),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
