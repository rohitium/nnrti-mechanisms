#!/usr/bin/env python3
"""Plot WT-vs-mutation sidechain-DOR distance traces as per-mutation figures."""
from __future__ import annotations

import argparse
import json
import re
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


_AA3 = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "Q": "GLN",
    "E": "GLU",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
}


def _mutation_sort_key(m: str) -> tuple[int, str]:
    if m == "WT":
        return (0, m)
    if "+" in m:
        return (2, m)
    return (1, m)


def _extract_components(mutation: str) -> list[dict[str, object]]:
    comps: list[dict[str, object]] = []
    for token in mutation.split("+"):
        tok = token.strip()
        m = re.match(r"^([A-Z])(\d+)([A-Z])$", tok)
        if not m:
            continue
        wt_aa, pos, mut_aa = m.groups()
        comps.append(
            {
                "token": tok,
                "wt_aa": wt_aa,
                "position": int(pos),
                "mut_aa": mut_aa,
                "wt_resname": _AA3.get(wt_aa),
                "mut_resname": _AA3.get(mut_aa),
            }
        )
    return comps


def _remap_to_local_workspace(candidate: Path | None, repo_root: Path) -> Path | None:
    if candidate is None:
        return None
    if candidate.exists():
        return candidate
    marker = "nnrti-mechanisms/"
    text = str(candidate)
    if marker not in text:
        return candidate
    rel = text.split(marker, 1)[1]
    mapped = repo_root / rel
    if mapped.exists():
        return mapped
    return candidate


def _replicate_inputs(row: pd.Series, repo_root: Path) -> tuple[Path, Path]:
    data = json.loads(Path(row["output_json"]).read_text())
    topo = Path(str(data.get("analysis_topology_pdb") or "").strip())
    dcd = Path(str(data.get("analysis_dcd") or "").strip())
    topo = _remap_to_local_workspace(topo, repo_root)
    dcd = _remap_to_local_workspace(dcd, repo_root)
    if topo is None or dcd is None or not topo.exists() or not dcd.exists():
        raise FileNotFoundError(f"Missing analysis files for {row['mutation']} rep{int(row['replicate'])}")
    return topo, dcd


def _choose_residue(universe, position: int, ligand_sel: str, expected_resname: str | None):
    from MDAnalysis.lib.distances import distance_array

    lig = universe.select_atoms(ligand_sel)
    if lig.n_atoms == 0:
        raise ValueError(f"DOR selection returned 0 atoms: {ligand_sel}")

    residues = list(universe.select_atoms(f"protein and resid {position} and not name H*").residues)
    if not residues:
        raise ValueError(f"No residue found for position {position}")

    if expected_resname:
        match = [r for r in residues if str(r.resname).upper() == str(expected_resname).upper()]
        if match:
            residues = match

    universe.trajectory[0]
    best = None
    best_d = np.inf
    for residue in residues:
        heavy = residue.atoms.select_atoms("not name H*")
        if heavy.n_atoms == 0:
            continue
        d = float(distance_array(heavy.positions, lig.positions, box=universe.dimensions).min())
        if d < best_d:
            best_d = d
            best = residue
    if best is None:
        raise ValueError(f"Could not select heavy-atom residue for position {position}")
    return best


def _sidechain_atoms(residue):
    sc = residue.atoms.select_atoms("not name N CA C O OXT and not name H*")
    if sc.n_atoms == 0:
        sc = residue.atoms.select_atoms("not name H*")
    return sc


def _interp_mean_trace(df: pd.DataFrame, x_col: str = "time_ns", y_col: str = "distance_angstrom", n_grid: int = 200):
    if df.empty:
        return None, None
    xmin = float(pd.to_numeric(df[x_col], errors="coerce").min())
    xmax = float(pd.to_numeric(df[x_col], errors="coerce").max())
    if not np.isfinite(xmin) or not np.isfinite(xmax) or xmax <= xmin:
        return None, None
    grid = np.linspace(xmin, xmax, n_grid)
    ys = []
    for _rep, grp in df.groupby("replicate"):
        g = grp.sort_values(x_col)
        x = g[x_col].to_numpy(dtype=float)
        y = g[y_col].to_numpy(dtype=float)
        if len(x) < 2:
            continue
        keep = np.r_[True, np.diff(x) > 0]
        x = x[keep]
        y = y[keep]
        if len(x) < 2:
            continue
        yi = np.interp(grid, x, y, left=np.nan, right=np.nan)
        yi[(grid < x.min()) | (grid > x.max())] = np.nan
        ys.append(yi)
    if not ys:
        return None, None
    return grid, np.nanmean(np.vstack(ys), axis=0)


def _metric_specs(comps: list[dict[str, object]]) -> list[tuple[str, str]]:
    if len(comps) == 1:
        return [("c1_to_dor", str(comps[0]["token"]))]
    c1 = str(comps[0]["token"])
    c2 = str(comps[1]["token"])
    return [
        ("c1_to_dor", f"{c1} to DOR"),
        ("c2_to_dor", f"{c2} to DOR"),
        ("c1_to_c2", f"{c1} to {c2}"),
    ]


def _trace_features(trace_df: pd.DataFrame) -> dict[str, float]:
    g = trace_df.sort_values("time_ns")
    x = pd.to_numeric(g["time_ns"], errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(g["distance_angstrom"], errors="coerce").to_numpy(dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    x = x[keep]
    y = y[keep]

    if x.size < 4:
        return {
            "n_points": int(x.size),
            "max_step_abs_angstrom": np.nan,
            "p95_step_abs_angstrom": np.nan,
            "late_max_step_abs_angstrom": np.nan,
            "p95_curvature_abs_angstrom": np.nan,
            "endpoint_shift_abs_angstrom": np.nan,
            "trace_range_angstrom": np.nan,
            "trace_std_angstrom": np.nan,
            "max_step_time_ns": np.nan,
        }

    dy = np.diff(y)
    dx = np.diff(x)
    valid = dx > 0
    dy = dy[valid]
    x_step = x[1:][valid]
    if dy.size == 0:
        return {
            "n_points": int(x.size),
            "max_step_abs_angstrom": np.nan,
            "p95_step_abs_angstrom": np.nan,
            "late_max_step_abs_angstrom": np.nan,
            "p95_curvature_abs_angstrom": np.nan,
            "endpoint_shift_abs_angstrom": np.nan,
            "trace_range_angstrom": float(np.nanmax(y) - np.nanmin(y)),
            "trace_std_angstrom": float(np.nanstd(y)),
            "max_step_time_ns": np.nan,
        }

    abs_step = np.abs(dy)
    max_idx = int(np.nanargmax(abs_step))
    xmax = float(np.nanmax(x))
    late_mask = x_step >= (0.7 * xmax)
    late_max = float(np.nanmax(abs_step[late_mask])) if np.any(late_mask) else float(np.nanmax(abs_step))

    d2 = np.diff(y, n=2)
    p95_curv = float(np.nanpercentile(np.abs(d2), 95)) if d2.size > 0 else np.nan

    n = len(y)
    edge_n = max(3, int(round(0.1 * n)))
    head_mean = float(np.nanmean(y[:edge_n]))
    tail_mean = float(np.nanmean(y[-edge_n:]))

    return {
        "n_points": int(x.size),
        "max_step_abs_angstrom": float(np.nanmax(abs_step)),
        "p95_step_abs_angstrom": float(np.nanpercentile(abs_step, 95)),
        "late_max_step_abs_angstrom": late_max,
        "p95_curvature_abs_angstrom": p95_curv,
        "endpoint_shift_abs_angstrom": float(abs(tail_mean - head_mean)),
        "trace_range_angstrom": float(np.nanmax(y) - np.nanmin(y)),
        "trace_std_angstrom": float(np.nanstd(y)),
        "max_step_time_ns": float(x_step[max_idx]),
    }


def _positive_robust_z(values: pd.Series) -> pd.Series:
    s = pd.to_numeric(values, errors="coerce")
    valid = s.dropna()
    out = pd.Series(np.zeros(len(s), dtype=float), index=s.index)
    if valid.empty:
        return out
    if valid.nunique() <= 1:
        return out

    med = float(valid.median())
    mad = float(np.median(np.abs(valid.to_numpy(dtype=float) - med)))
    if mad > 1e-12:
        z = (s - med) / (1.4826 * mad)
    else:
        std = float(valid.std(ddof=0))
        if std <= 1e-12:
            return out
        z = (s - float(valid.mean())) / std
    out = z.clip(lower=0.0, upper=8.0).fillna(0.0)
    return out


def curate_interesting_traces(
    timeseries_df: pd.DataFrame,
    output_csv: Path,
    plots_dir: Path,
    top_n: int = 100,
    min_score: float = 5.0,
) -> pd.DataFrame:
    group_cols = ["mutation", "system", "metric", "replicate"]
    rows: list[dict[str, object]] = []

    for keys, grp in timeseries_df.groupby(group_cols, dropna=False):
        mutation, system, metric, replicate = keys
        feat = _trace_features(grp)
        first = grp.iloc[0]
        rows.append(
            {
                "mutation": str(mutation),
                "system": str(system),
                "metric": str(metric),
                "replicate": int(replicate),
                "safe_label": str(first.get("safe_label", "")),
                "output_json": str(first.get("output_json", "")),
                "analysis_dcd": str(first.get("analysis_dcd", "")),
                "analysis_topology_pdb": str(first.get("analysis_topology_pdb", "")),
                "drm_plot_path": str(
                    plots_dir / f"{str(mutation).replace('+', '_')}_drm_distance_timeseries.png"
                ),
                **feat,
            }
        )

    trace_df = pd.DataFrame(rows)
    if trace_df.empty:
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        trace_df.to_csv(output_csv, index=False)
        return trace_df

    panel_cols = ["mutation", "system", "metric"]
    z_specs = [
        ("max_step_abs_angstrom", "z_max_step"),
        ("late_max_step_abs_angstrom", "z_late_jump"),
        ("p95_curvature_abs_angstrom", "z_curvature"),
        ("endpoint_shift_abs_angstrom", "z_endpoint_shift"),
    ]
    for src, dest in z_specs:
        trace_df[dest] = (
            trace_df.groupby(panel_cols, dropna=False)[src]
            .transform(_positive_robust_z)
            .astype(float)
        )

    trace_df["interesting_score"] = (
        2.0 * trace_df["z_max_step"]
        + 1.6 * trace_df["z_late_jump"]
        + 1.3 * trace_df["z_curvature"]
        + 1.0 * trace_df["z_endpoint_shift"]
    )

    reasons: list[str] = []
    reason_cols: list[list[str]] = []
    for _, row in trace_df.iterrows():
        r: list[str] = []
        if float(row.get("max_step_abs_angstrom", np.nan)) >= 2.6:
            r.append("hard_jump>=2.6A")
        if float(row.get("late_max_step_abs_angstrom", np.nan)) >= 2.4:
            r.append("late_jump>=2.4A")
        if float(row.get("p95_curvature_abs_angstrom", np.nan)) >= 2.4:
            r.append("bumpy_curvature>=2.4A")
        if float(row.get("endpoint_shift_abs_angstrom", np.nan)) >= 1.0:
            r.append("endpoint_shift>=1.0A")
        if float(row.get("z_max_step", 0.0)) >= 2.5:
            r.append("panel_outlier_step")
        if float(row.get("z_curvature", 0.0)) >= 2.5:
            r.append("panel_outlier_bumpy")
        reason_cols.append(r)
        reasons.append(";".join(r))
    trace_df["interesting_reasons"] = reasons
    trace_df["n_reasons"] = [len(r) for r in reason_cols]

    has_dynamic_signal = (
        (trace_df["z_max_step"] >= 1.5)
        | (trace_df["z_late_jump"] >= 1.5)
        | (trace_df["z_curvature"] >= 1.5)
        | (trace_df["endpoint_shift_abs_angstrom"] >= 1.0)
    )
    trace_df["is_interesting"] = (
        ((trace_df["interesting_score"] >= float(min_score)) & has_dynamic_signal)
        | (trace_df["n_reasons"] > 0)
    )
    trace_df = trace_df.sort_values(
        ["is_interesting", "interesting_score", "n_reasons", "max_step_abs_angstrom"],
        ascending=[False, False, False, False],
    ).reset_index(drop=True)
    if top_n is not None and int(top_n) > 0 and len(trace_df) > int(top_n):
        trace_df = trace_df.iloc[: int(top_n)].copy()

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    trace_df.to_csv(output_csv, index=False)
    return trace_df


def _infer_total_ns_from_output_json(output_json_path: Path) -> float | None:
    state_csv = None
    m = re.match(r"^(.+)_rep(\d{2})\.json$", output_json_path.name)
    if m:
        safe = m.group(1)
        rep = int(m.group(2))
        state_csv = output_json_path.parent / f"{safe}_rep{rep:02d}_md_state.csv"
    if state_csv is None or not state_csv.exists():
        return None
    try:
        sdf = pd.read_csv(state_csv)
    except Exception:
        return None
    step_col = None
    for c in ('#"Step"', "Step"):
        if c in sdf.columns:
            step_col = c
            break
    if step_col is None or sdf.empty:
        return None
    steps = pd.to_numeric(sdf[step_col], errors="coerce").dropna()
    if steps.empty:
        return None
    return float(steps.max()) * 2.0 / 1_000_000.0


def _collect_system_rows(
    rows_df: pd.DataFrame,
    system_label: str,
    comps: list[dict[str, object]],
    ligand_resname: str,
    frame_stride: int,
    repo_root: Path,
    resid_offset: int,
) -> list[dict]:
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import distance_array

    out: list[dict] = []
    for _, row in rows_df.sort_values("replicate").iterrows():
        replicate = int(row["replicate"])
        try:
            topo, dcd = _replicate_inputs(row, repo_root)
        except FileNotFoundError:
            continue  # replicate not yet complete; skip silently
        u = mda.Universe(str(topo), str(dcd))
        dor = u.select_atoms(f"resname {ligand_resname} and not name H*")
        if dor.n_atoms == 0:
            continue

        if system_label == "WT":
            c1_resname = str(comps[0]["wt_resname"] or "")
            c2_resname = str(comps[1]["wt_resname"] or "") if len(comps) > 1 else None
        else:
            c1_resname = str(comps[0]["mut_resname"] or "")
            c2_resname = str(comps[1]["mut_resname"] or "") if len(comps) > 1 else None

        p1 = int(comps[0]["position"]) + int(resid_offset)
        r1 = _choose_residue(
            u, p1, f"resname {ligand_resname} and not name H*", c1_resname or None
        )
        sc1 = _sidechain_atoms(r1)
        if sc1.n_atoms == 0:
            continue

        sc2 = None
        if len(comps) > 1:
            p2 = int(comps[1]["position"]) + int(resid_offset)
            r2 = _choose_residue(
                u, p2, f"resname {ligand_resname} and not name H*", c2_resname or None
            )
            sc2 = _sidechain_atoms(r2)
            if sc2.n_atoms == 0:
                continue

        max_frame = max(1, len(u.trajectory) - 1)
        total_ns = _infer_total_ns_from_output_json(Path(str(row["output_json"])))
        if total_ns is None or not np.isfinite(total_ns) or total_ns <= 0:
            total_ns = 2.0
        for ts in u.trajectory[:: max(1, frame_stride)]:
            t_ns = (float(ts.frame) / float(max_frame)) * float(total_ns)
            d1 = float(distance_array(sc1.positions, dor.positions, box=u.dimensions).min())
            out.append(
                {
                    "mutation": str(row["mutation"]),
                    "safe_label": str(row.get("safe_label", "")),
                    "system": system_label,
                    "replicate": replicate,
                    "time_ns": t_ns,
                    "metric": "c1_to_dor",
                    "distance_angstrom": d1,
                    "output_json": str(row.get("output_json", "")),
                    "analysis_dcd": str(dcd),
                    "analysis_topology_pdb": str(topo),
                }
            )
            if sc2 is not None:
                d2 = float(distance_array(sc2.positions, dor.positions, box=u.dimensions).min())
                d12 = float(distance_array(sc1.positions, sc2.positions, box=u.dimensions).min())
                out.append(
                    {
                        "mutation": str(row["mutation"]),
                        "safe_label": str(row.get("safe_label", "")),
                        "system": system_label,
                        "replicate": replicate,
                        "time_ns": t_ns,
                        "metric": "c2_to_dor",
                        "distance_angstrom": d2,
                        "output_json": str(row.get("output_json", "")),
                        "analysis_dcd": str(dcd),
                        "analysis_topology_pdb": str(topo),
                    }
                )
                out.append(
                    {
                        "mutation": str(row["mutation"]),
                        "safe_label": str(row.get("safe_label", "")),
                        "system": system_label,
                        "replicate": replicate,
                        "time_ns": t_ns,
                        "metric": "c1_to_c2",
                        "distance_angstrom": d12,
                        "output_json": str(row.get("output_json", "")),
                        "analysis_dcd": str(dcd),
                        "analysis_topology_pdb": str(topo),
                    }
                )
    return out


def main() -> int:
    warnings.filterwarnings(
        "ignore",
        category=DeprecationWarning,
        message=r"DCDReader currently makes independent timesteps.*",
    )
    parser = argparse.ArgumentParser(description="Plot sidechain-DOR distances for all mutations vs WT.")
    parser.add_argument("--manifest", type=Path, default=Path("manifests/md_manifest.csv"))
    parser.add_argument("--ligand-resname", type=str, default="2KW")
    parser.add_argument("--frame-stride", type=int, default=1)
    parser.add_argument("--resid-offset", type=int, default=-3)
    parser.add_argument("--plots-dir", type=Path, default=Path("results/plots/drm_distances"))
    parser.add_argument("--output-csv", type=Path, default=Path("results/drm_sidechain_distance_timeseries_all_mutations.csv"))
    parser.add_argument(
        "--interesting-csv",
        type=Path,
        default=Path("results/drm_sidechain_distance_interesting_traces.csv"),
    )
    parser.add_argument("--interesting-top-n", type=int, default=60)
    parser.add_argument("--interesting-min-score", type=float, default=5.0)
    args = parser.parse_args()

    if not args.manifest.exists():
        raise FileNotFoundError(args.manifest)

    repo_root = Path(__file__).resolve().parents[3]
    mf = pd.read_csv(args.manifest)
    wt_df = mf[mf["mutation"] == "WT"].copy()
    mut_df = mf[mf["mutation"] != "WT"].copy()
    if wt_df.empty or mut_df.empty:
        raise ValueError("Manifest must contain WT and non-WT mutations.")

    muts = sorted(mut_df["mutation"].unique(), key=lambda m: _mutation_sort_key(str(m)))

    # Pre-collect all rows so we can size the panel stack.
    by_mut_data: dict[str, pd.DataFrame] = {}
    metric_titles: dict[str, list[tuple[str, str]]] = {}
    all_rows: list[dict] = []

    for mut in muts:
        comps = _extract_components(str(mut))
        if not comps:
            continue
        mdf = mut_df[mut_df["mutation"] == mut].copy()
        wt_lookup = wt_df[wt_df["replicate"].isin(mdf["replicate"])].copy()
        if wt_lookup.empty:
            wt_lookup = wt_df.copy()

        rows = []
        rows.extend(
            _collect_system_rows(
                mdf, "Mutant", comps, args.ligand_resname, args.frame_stride, repo_root, args.resid_offset
            )
        )
        rows.extend(
            _collect_system_rows(
                wt_lookup, "WT", comps, args.ligand_resname, args.frame_stride, repo_root, args.resid_offset
            )
        )
        if not rows:
            continue
        block = pd.DataFrame(rows)
        block["mutation"] = mut
        by_mut_data[str(mut)] = block
        metric_titles[str(mut)] = _metric_specs(comps)
        all_rows.extend(block.to_dict(orient="records"))

    if not by_mut_data:
        raise ValueError("No sidechain-DOR traces were generated.")

    mut_color = "#1f77b4"
    wt_color = "#444444"
    args.plots_dir.mkdir(parents=True, exist_ok=True)

    for mut in muts:
        if mut not in by_mut_data:
            continue
        block = by_mut_data[mut]
        specs = metric_titles[mut]
        nrows = len(specs)
        fig, axes = plt.subplots(nrows, 1, figsize=(9.0, 3.2 * nrows), squeeze=False)
        axes_list = axes[:, 0].tolist()

        for i, (metric_key, panel_title) in enumerate(specs):
            ax = axes_list[i]
            sub = block[block["metric"] == metric_key].copy()
            if sub.empty:
                ax.set_visible(False)
                continue

            # Light replicate trajectories
            for system, color in [("WT", wt_color), ("Mutant", mut_color)]:
                ss = sub[sub["system"] == system]
                for _rep, grp in ss.groupby("replicate"):
                    g = grp.sort_values("time_ns")
                    ax.plot(
                        g["time_ns"].to_numpy(dtype=float),
                        g["distance_angstrom"].to_numpy(dtype=float),
                        color=color,
                        alpha=0.25,
                        linewidth=0.8,
                    )

            # Mean trace + dashed global mean for each system
            for system, color in [("WT", wt_color), ("Mutant", mut_color)]:
                ss = sub[sub["system"] == system]
                x_mean, y_mean = _interp_mean_trace(ss)
                if x_mean is not None:
                    ax.plot(
                        x_mean,
                        y_mean,
                        color=color,
                        linewidth=2.0,
                        alpha=0.95,
                        label=system,
                    )
                if not ss.empty and ss["distance_angstrom"].notna().any():
                    mean_val = float(ss["distance_angstrom"].mean())
                    ax.axhline(mean_val, color=color, linestyle="--", linewidth=1.0, alpha=0.9)

            # Titles: single DRM uses mutation name only.
            if len(specs) == 1:
                ax.set_title(str(mut), fontsize=10, fontweight="bold")
            else:
                ax.set_title(panel_title.replace("ligand", "DOR"), fontsize=9)

            xmax = float(pd.to_numeric(sub["time_ns"], errors="coerce").max())
            if np.isfinite(xmax) and xmax > 0:
                ax.set_xlim(0.0, xmax)
            ax.set_xlabel("Time (ns)")
            ax.set_ylabel("Min distance (Å)")
            ax.grid(alpha=0.25, linestyle=":")
            ax.legend(frameon=False, fontsize=8, loc="best")

        fig.suptitle(f"{mut}: Sidechain-DOR Distance Plots", fontsize=11, fontweight="bold", y=0.995)
        fig.tight_layout()
        out_plot = args.plots_dir / f"{str(mut).replace('+', '_')}_drm_distance_timeseries.png"
        fig.savefig(out_plot, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Wrote {out_plot}")

    out_df = pd.DataFrame(all_rows)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.output_csv, index=False)
    print(f"Wrote {args.output_csv}")
    interesting_df = curate_interesting_traces(
        out_df,
        output_csv=args.interesting_csv,
        plots_dir=args.plots_dir,
        top_n=int(args.interesting_top_n),
        min_score=float(args.interesting_min_score),
    )
    print(f"Wrote {args.interesting_csv} (rows={len(interesting_df)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
