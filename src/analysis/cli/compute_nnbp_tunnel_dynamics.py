#!/usr/bin/env python3
"""Compute NNBP tunnel gate distances over MD trajectories.

Measures two proxy distances that characterize the entrance/aperture of the
Non-Nucleoside Inhibitor Binding Pocket (NNBP) tunnel:

  gate1 (vertical):  K101(Cα) ↔ Y188(Cα)  — upper-wall to floor
  gate2 (lateral):   V106(Cβ) ↔ Y181(Cβ)  — lateral walls

Works for both holo (DOR-bound) and apo (ligand-free) simulations, making it
the key observable for Hypothesis 2 (kinetic tunnel-opening mechanism).

Output CSVs: results/nnbp_tunnel_dynamics.csv  (per-frame rows)
             results/nnbp_tunnel_summary.csv    (per-replicate mean ± std)
"""
from __future__ import annotations

import argparse
import json
import logging
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd


_NNBP_GATE_RESIDUES = {
    # (resname_wt, position_4ncg_numbering)
    "K101": ("LYS", 101),
    "Y188": ("TYR", 188),
    "V106": ("VAL", 106),
    "Y181": ("TYR", 181),
    # Additional floor/ceiling markers
    "F227": ("PHE", 227),
    "P225": ("PRO", 225),
    "L100": ("LEU", 100),
    "K103": ("LYS", 103),
}

# Gate pair definitions: (label, res_a_key, atom_a, res_b_key, atom_b)
_GATE_PAIRS = [
    ("gate_K101_Y188_CA",  "K101", "CA",  "Y188", "CA"),
    ("gate_V106_Y181_CB",  "V106", "CB",  "Y181", "CB"),
    ("gate_L100_F227_CA",  "L100", "CA",  "F227", "CA"),
    ("gate_K103_P225_CA",  "K103", "CA",  "P225", "CA"),
]


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
        raise FileNotFoundError(
            f"Missing analysis files for {row['mutation']} rep{int(row['replicate'])}"
        )
    return topo, dcd


def _infer_total_ns(output_json_path: Path) -> float:
    """Read production length from JSON or state CSV; default 100 ns."""
    # Try JSON key directly
    try:
        j = json.loads(output_json_path.read_text())
        steps = int(j.get("md_production_steps_completed") or j.get("md_production_steps") or 0)
        if steps > 0:
            return steps * 2.0 / 1_000_000.0
    except Exception:
        pass
    # Fall back to state CSV
    m = re.match(r"^(.+)_rep(\d{2})\.json$", output_json_path.name)
    if m:
        state_csv = output_json_path.parent / f"{m.group(1)}_rep{m.group(2)}_md_state.csv"
        if state_csv.exists():
            try:
                sdf = pd.read_csv(state_csv)
                for col in ('#"Step"', "Step"):
                    if col in sdf.columns:
                        steps = pd.to_numeric(sdf[col], errors="coerce").dropna()
                        if not steps.empty:
                            return float(steps.max()) * 2.0 / 1_000_000.0
            except Exception:
                pass
    return 100.0


def _select_atom(universe, position: int, resid_offset: int, atom_name: str):
    """Select a single named atom at a given sequence position."""
    import MDAnalysis as mda  # noqa: F401

    resid = position + resid_offset
    sel = universe.select_atoms(f"protein and resid {resid} and name {atom_name}")
    if sel.n_atoms == 0:
        # Some residues (Pro, Gly) lack CB — fall back to CA
        if atom_name == "CB":
            sel = universe.select_atoms(f"protein and resid {resid} and name CA")
    if sel.n_atoms == 0:
        raise ValueError(f"No atom {atom_name} at resid {resid} (position {position})")
    # Return just the first atom if multiple matches (shouldn't happen in stripped topology)
    return sel[0:1]


def _process_replicate(
    row: pd.Series,
    repo_root: Path,
    resid_offset: int,
    frame_stride: int,
) -> list[dict]:
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import distance_array

    topo, dcd = _replicate_inputs(row, repo_root)
    u = mda.Universe(str(topo), str(dcd))
    mutation = str(row["mutation"])
    replicate = int(row["replicate"])
    total_ns = _infer_total_ns(Path(str(row["output_json"])))
    n_frames = len(u.trajectory)

    # Pre-select atoms for all gate pairs
    atom_pairs: list[tuple[str, object, object]] = []
    for label, key_a, atom_a, key_b, atom_b in _GATE_PAIRS:
        pos_a = _NNBP_GATE_RESIDUES[key_a][1]
        pos_b = _NNBP_GATE_RESIDUES[key_b][1]
        try:
            ag_a = _select_atom(u, pos_a, resid_offset, atom_a)
            ag_b = _select_atom(u, pos_b, resid_offset, atom_b)
        except ValueError as exc:
            logging.warning(f"  {mutation} rep{replicate}: skipping {label} — {exc}")
            continue
        atom_pairs.append((label, ag_a, ag_b))

    if not atom_pairs:
        return []

    out: list[dict] = []
    for ts in u.trajectory[:: max(1, frame_stride)]:
        time_ns = (float(ts.frame) / max(1, n_frames - 1)) * total_ns
        row_base = {
            "mutation": mutation,
            "safe_label": str(row.get("safe_label", "")),
            "replicate": replicate,
            "frame": int(ts.frame),
            "time_ns": time_ns,
        }
        for label, ag_a, ag_b in atom_pairs:
            d = float(
                distance_array(ag_a.positions, ag_b.positions, box=u.dimensions).min()
            )
            out.append({**row_base, "gate": label, "distance_angstrom": d})

    return out


def main() -> int:
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    parser = argparse.ArgumentParser(
        description="Compute NNBP tunnel gate distances over trajectories."
    )
    parser.add_argument("--manifest", type=Path, default=Path("results/md_manifest.csv"))
    parser.add_argument("--resid-offset", type=int, default=-3)
    parser.add_argument("--frame-stride", type=int, default=1)
    parser.add_argument(
        "--mutations",
        nargs="*",
        default=None,
        help="Subset of mutations to process (default: all in manifest)",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=Path("results/nnbp_tunnel_dynamics.csv"),
    )
    parser.add_argument(
        "--summary-csv",
        type=Path,
        default=Path("results/nnbp_tunnel_summary.csv"),
    )
    parser.add_argument("--plots-dir", type=Path, default=Path("results/plots/nnbp_tunnel"))
    args = parser.parse_args()

    if not args.manifest.exists():
        raise FileNotFoundError(args.manifest)

    repo_root = Path(__file__).resolve().parents[3]
    mf = pd.read_csv(args.manifest)

    if args.mutations:
        mf = mf[mf["mutation"].isin(args.mutations)].copy()
        if mf.empty:
            logging.error(f"No manifest rows match --mutations {args.mutations}")
            return 1

    all_rows: list[dict] = []
    for _, row in mf.iterrows():
        mut = str(row["mutation"])
        rep = int(row["replicate"])
        logging.info(f"Processing {mut} rep{rep}...")
        try:
            rows = _process_replicate(row, repo_root, args.resid_offset, args.frame_stride)
            all_rows.extend(rows)
            logging.info(f"  → {len(rows)} frame-gate rows")
        except Exception as exc:
            logging.error(f"  FAILED {mut} rep{rep}: {exc}")

    if not all_rows:
        logging.error("No tunnel dynamics data collected.")
        return 1

    out_df = pd.DataFrame(all_rows)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.output_csv, index=False)
    logging.info(f"Wrote {args.output_csv} ({len(out_df)} rows)")

    # Summary: per-mutation per-gate mean ± std
    summary = (
        out_df.groupby(["mutation", "gate"])["distance_angstrom"]
        .agg(mean="mean", std="std", n="count")
        .reset_index()
    )
    args.summary_csv.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.summary_csv, index=False)
    logging.info(f"Wrote {args.summary_csv}")

    # Plots: per-gate, all mutations overlaid
    try:
        _plot_tunnel_dynamics(out_df, args.plots_dir)
    except Exception as exc:
        logging.error(f"Plotting failed: {exc}")

    return 0


def _plot_tunnel_dynamics(df: pd.DataFrame, plots_dir: Path) -> None:
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm

    plots_dir.mkdir(parents=True, exist_ok=True)
    gates = df["gate"].unique()

    for gate in gates:
        sub = df[df["gate"] == gate].copy()
        mutations = sorted(sub["mutation"].unique())
        colors = cm.tab20(np.linspace(0, 1, len(mutations)))

        fig, ax = plt.subplots(figsize=(10, 4))
        for mut, color in zip(mutations, colors):
            ms = sub[sub["mutation"] == mut]
            for _, grp in ms.groupby("replicate"):
                g = grp.sort_values("time_ns")
                ax.plot(
                    g["time_ns"].to_numpy(float),
                    g["distance_angstrom"].to_numpy(float),
                    color=color,
                    alpha=0.4,
                    linewidth=0.8,
                )
            # Mean trace per mutation
            by_time = ms.groupby("time_ns")["distance_angstrom"].mean()
            ax.plot(
                by_time.index.to_numpy(float),
                by_time.to_numpy(float),
                color=color,
                linewidth=1.5,
                label=mut,
            )

        gate_label = gate.replace("gate_", "").replace("_", " ")
        ax.set_title(f"NNBP Tunnel Gate: {gate_label}", fontsize=11, fontweight="bold")
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("Distance (Å)")
        ax.grid(alpha=0.25, linestyle=":")
        ax.legend(frameon=False, fontsize=7, bbox_to_anchor=(1.01, 1), loc="upper left")
        fig.tight_layout()

        out_path = plots_dir / f"{gate}_timeseries.png"
        fig.savefig(out_path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        logging.info(f"Wrote {out_path}")

    # Distribution plot: compare WT vs resistant mutants for primary gate
    primary = "gate_K101_Y188_CA"
    if primary in gates:
        _plot_gate_distributions(df[df["gate"] == primary], primary, plots_dir)


def _plot_gate_distributions(df: pd.DataFrame, gate: str, plots_dir: Path) -> None:
    import matplotlib.pyplot as plt

    mutations = sorted(df["mutation"].unique(), key=lambda m: (m != "WT", "+" in m, m))
    fig, ax = plt.subplots(figsize=(max(6, len(mutations) * 0.6), 4))

    data = [
        df[df["mutation"] == mut]["distance_angstrom"].dropna().to_numpy(float)
        for mut in mutations
    ]
    positions = list(range(1, len(mutations) + 1))
    bp = ax.boxplot(data, positions=positions, patch_artist=True, widths=0.5)

    wt_idx = mutations.index("WT") if "WT" in mutations else -1
    for i, patch in enumerate(bp["boxes"]):
        patch.set_facecolor("#2196F3" if i == wt_idx else "#FF5722")
        patch.set_alpha(0.7)

    ax.set_xticks(positions)
    ax.set_xticklabels(mutations, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Gate distance K101–Y188 Cα (Å)")
    ax.set_title("NNBP Tunnel Aperture Distribution by Mutation", fontweight="bold")
    ax.grid(axis="y", alpha=0.3, linestyle=":")
    fig.tight_layout()

    out_path = plots_dir / f"{gate}_distribution_by_mutation.png"
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    logging.info(f"Wrote {out_path}")


if __name__ == "__main__":
    raise SystemExit(main())
