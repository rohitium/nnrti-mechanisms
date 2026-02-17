#!/usr/bin/env python3
"""Plot NNRTI pocket-volume distributions by mutation.

Inputs:
  - results/pocket_volume_profiles.csv (frame-level pocket volume metric)
  - results/md_manifest.csv (fold reductions, mutation ordering)

Outputs:
  - results/plots/pocket_volume_distribution_by_mutation.png (default)
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def _mutation_order(manifest_csv: Path) -> list[str]:
    m = pd.read_csv(manifest_csv)
    m["fold_reduction"] = pd.to_numeric(m.get("fold_reduction"), errors="coerce")
    uniq = (
        m[["mutation", "fold_reduction"]]
        .drop_duplicates()
        .sort_values(["fold_reduction", "mutation"], na_position="first")
    )
    muts = uniq["mutation"].astype(str).tolist()
    # Prefer WT first.
    if "WT" in muts:
        muts = ["WT"] + [x for x in muts if x != "WT"]
    return muts


def plot_pocket_volume_distributions(
    pocket_volume_profiles_csv: Path,
    manifest_csv: Path,
    md_runs_dir: Path,
    output_png: Path,
    max_points_per_mutation: int,
    last_window_ns: float,
) -> None:
    import matplotlib.pyplot as plt

    df = pd.read_csv(pocket_volume_profiles_csv)
    if df.empty:
        raise ValueError(f"Empty pocket volume profiles: {pocket_volume_profiles_csv}")

    df["pocket_volume_proxy_angstrom3"] = pd.to_numeric(df["pocket_volume_proxy_angstrom3"], errors="coerce")
    df = df.dropna(subset=["mutation", "pocket_volume_proxy_angstrom3"]).copy()

    # Build a robust time axis from md_state.csv (preferred), then filter to the last window.
    def _infer_total_ns_from_state_csv(safe_label: str, replicate: int) -> float | None:
        state_csv = (
            md_runs_dir
            / str(safe_label)
            / f"rep_{int(replicate):02d}"
            / f"{str(safe_label)}_rep{int(replicate):02d}_md_state.csv"
        )
        if not state_csv.exists():
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
        max_step = pd.to_numeric(sdf[step_col], errors="coerce").dropna()
        if max_step.empty:
            return None
        return float(max_step.max()) * 2.0 / 1_000_000.0  # 2 fs timestep

    if {"safe_label", "replicate", "frame_index"}.issubset(df.columns):
        parts = []
        for (safe, rep), sub in df.groupby(["safe_label", "replicate"], dropna=False):
            sub = sub.copy()
            total_ns = _infer_total_ns_from_state_csv(str(safe), int(rep))
            frame = pd.to_numeric(sub["frame_index"], errors="coerce")
            max_frame = float(frame.max()) if frame.notna().any() else np.nan
            if total_ns is not None and np.isfinite(total_ns) and total_ns > 0 and np.isfinite(max_frame) and max_frame > 0:
                sub["time_ns"] = (frame / max_frame) * float(total_ns)
            else:
                sub["time_ns"] = pd.to_numeric(sub.get("time_ps", np.nan), errors="coerce") / 1000.0
            parts.append(sub)
        df = pd.concat(parts, ignore_index=True) if parts else df
    else:
        df["time_ns"] = pd.to_numeric(df.get("time_ps", np.nan), errors="coerce") / 1000.0

    if float(last_window_ns) > 0.0 and "time_ns" in df.columns and {"safe_label", "replicate"}.issubset(df.columns):
        max_t = df.groupby(["safe_label", "replicate"], dropna=False)["time_ns"].transform("max")
        can_filter = df["time_ns"].notna() & max_t.notna()
        df = df[~can_filter | (df["time_ns"] >= (max_t - float(last_window_ns)))].copy()

    order = _mutation_order(manifest_csv)
    order = [m for m in order if m in set(df["mutation"].astype(str))]
    if not order:
        raise ValueError("No overlapping mutations between manifest and pocket_volume_profiles.csv.")

    # Downsample for scatter overlay (keeps plot readable).
    rng = np.random.default_rng(0)
    parts = []
    for mut in order:
        g = df[df["mutation"].astype(str) == mut]
        if g.empty:
            continue
        if len(g) > max_points_per_mutation:
            take = rng.choice(g.index.to_numpy(), size=max_points_per_mutation, replace=False)
            parts.append(g.loc[take])
        else:
            parts.append(g)
    scatter_df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()

    data = [df[df["mutation"].astype(str) == mut]["pocket_volume_proxy_angstrom3"].to_numpy(dtype=float) for mut in order]

    fig, ax = plt.subplots(figsize=(14, 5))
    bp = ax.boxplot(
        data,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 1.3},
        whiskerprops={"linewidth": 1.0},
        capprops={"linewidth": 1.0},
    )
    for patch in bp["boxes"]:
        patch.set_facecolor("#cfe8f3")
        patch.set_edgecolor("#4d4d4d")
        patch.set_alpha(0.9)

    if not scatter_df.empty:
        x_positions = {mut: i + 1 for i, mut in enumerate(order)}
        xs = []
        ys = []
        for mut, g in scatter_df.groupby(scatter_df["mutation"].astype(str)):
            x0 = x_positions.get(mut)
            if x0 is None:
                continue
            jitter = rng.normal(0.0, 0.06, size=len(g))
            xs.extend((x0 + jitter).tolist())
            ys.extend(g["pocket_volume_proxy_angstrom3"].to_numpy(dtype=float).tolist())
        ax.scatter(xs, ys, s=6, color="#1f77b4", alpha=0.25, linewidths=0)

    ax.set_xticks(range(1, len(order) + 1))
    ax.set_xticklabels(order, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("Pocket volume (Å$^3$)")
    title_suffix = f" (last {float(last_window_ns):g} ns)" if float(last_window_ns) > 0 else ""
    ax.set_title(f"NNRTI Pocket Volume Distributions Across Doravirine Resistance Panel{title_suffix}", fontsize=12, fontweight="bold")
    ax.grid(True, axis="y", alpha=0.25)

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output_png}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot pocket volume distributions by mutation.")
    parser.add_argument("--profiles", type=Path, default=Path("results/pocket_volume_profiles.csv"))
    parser.add_argument("--manifest", type=Path, default=Path("results/md_manifest.csv"))
    parser.add_argument("--md-runs-dir", type=Path, default=Path("results/md_runs"))
    parser.add_argument("--output", type=Path, default=Path("results/plots/pocket_volume_distribution_by_mutation.png"))
    parser.add_argument("--max-points", type=int, default=1200, help="Downsample scatter overlay per mutation.")
    parser.add_argument("--last-window-ns", type=float, default=1.0, help="Restrict to the last N ns of each trajectory.")
    args = parser.parse_args()

    if not args.profiles.exists():
        raise FileNotFoundError(f"Missing: {args.profiles}")
    if not args.manifest.exists():
        raise FileNotFoundError(f"Missing: {args.manifest}")

    plot_pocket_volume_distributions(
        pocket_volume_profiles_csv=args.profiles,
        manifest_csv=args.manifest,
        md_runs_dir=args.md_runs_dir,
        output_png=args.output,
        max_points_per_mutation=max(200, int(args.max_points)),
        last_window_ns=float(args.last_window_ns),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
