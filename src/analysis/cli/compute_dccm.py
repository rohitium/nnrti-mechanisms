#!/usr/bin/env python3
"""Compute Dynamic Cross-Correlation Matrix (DCCM) for NNBP vs polymerase domains.

For each replicate, computes the per-residue Cα displacement cross-correlation
matrix. Reports:
  1. Full NxN DCCM (saved as .npy)
  2. Mean cross-correlation between the NNBP residue set and the polymerase
     fingers/palm sub-domain, as a scalar per replicate (the allosteric coupling score)
  3. Heatmap plots

This directly tests Hypothesis 1 (F227C alone disrupts RT processivity via
coupling to fingers-palm domain movements) and Hypothesis 2 (epistatic partners
V106A/A98G restore NNBP–fingers coupling).

Usage:
    python -m src.analysis.cli.compute_dccm --manifest manifests/md_manifest.csv

Output:
    results/dccm/<mutation>_rep<N>_dccm.npy
    results/dccm_allosteric_coupling.csv  (per-replicate scalar coupling scores)
    results/plots/dccm/                   (heatmap PNGs)
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


# ── Residue range definitions (4NCG p66 numbering) ────────────────────────
# NNBP: residues directly lining the binding pocket
_NNBP_RESIDUES = list(range(98, 115)) + list(range(178, 192)) + list(range(223, 230))

# Polymerase fingers domain (β1–β6, roughly residues 1–85)
_FINGERS_RESIDUES = list(range(1, 86))

# Polymerase palm domain (β7–β10, connection loops, roughly 86–215 minus NNBP)
_PALM_RESIDUES = [r for r in range(86, 216) if r not in _NNBP_RESIDUES]

# Thumb domain (roughly 244–310 in p66)
_THUMB_RESIDUES = list(range(244, 311))


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


def _compute_dccm(positions: np.ndarray) -> np.ndarray:
    """Compute normalised cross-correlation matrix from (n_frames, n_atoms, 3) array.

    Returns NxN matrix with values in [-1, 1].
    """
    n_frames, n_atoms, _ = positions.shape
    # Mean-subtract per atom
    disp = positions - positions.mean(axis=0, keepdims=True)  # (F, N, 3)
    # Flatten xyz into per-atom displacement vectors
    disp_flat = disp.reshape(n_frames, n_atoms * 3)  # (F, N*3)

    # Per-atom magnitude array: sqrt( <Δr·Δr> )
    # We need <Δr_i · Δr_j> / sqrt(<|Δr_i|²> <|Δr_j|²>)
    # Compute per-atom squared magnitudes
    mag2 = np.sum(disp ** 2, axis=2)  # (F, N)
    mean_mag2 = mag2.mean(axis=0)     # (N,)

    # Cross-dot-product: <Δr_i · Δr_j>
    # Use matrix multiply over frames: disp (F,N,3), disp^T -> (N,N)
    # <Δr_i · Δr_j> = mean over frames of sum_k disp[f,i,k]*disp[f,j,k]
    cross = np.einsum("fid,fjd->ij", disp, disp) / n_frames  # (N, N)

    denom = np.sqrt(np.outer(mean_mag2, mean_mag2))
    denom = np.where(denom < 1e-12, 1.0, denom)
    dccm = cross / denom
    return np.clip(dccm, -1.0, 1.0)


def _mean_coupling(dccm: np.ndarray, idx_a: list[int], idx_b: list[int]) -> float:
    """Mean absolute cross-correlation between two sets of residue indices."""
    if not idx_a or not idx_b:
        return float("nan")
    sub = np.abs(dccm[np.ix_(idx_a, idx_b)])
    return float(sub.mean())


def _process_replicate(
    row: pd.Series,
    repo_root: Path,
    resid_offset: int,
    frame_stride: int,
    dccm_dir: Path,
) -> dict:
    import MDAnalysis as mda

    topo, dcd = _replicate_inputs(row, repo_root)
    mutation = str(row["mutation"])
    replicate = int(row["replicate"])
    safe = str(row.get("safe_label", mutation.replace("+", "_")))

    u = mda.Universe(str(topo), str(dcd))

    # Build protein Cα selection ordered by resid
    ca = u.select_atoms("protein and name CA")
    if ca.n_atoms == 0:
        raise ValueError("No Cα atoms found")

    resids = ca.resids  # shape (N,)

    # Collect Cα positions over trajectory
    pos_list: list[np.ndarray] = []
    for ts in u.trajectory[:: max(1, frame_stride)]:
        pos_list.append(ca.positions.copy())  # (N, 3)

    if len(pos_list) < 10:
        raise ValueError(f"Too few frames ({len(pos_list)}) for DCCM")

    positions = np.stack(pos_list, axis=0)  # (F, N, 3)
    dccm = _compute_dccm(positions)         # (N, N)

    # Save raw DCCM
    dccm_dir.mkdir(parents=True, exist_ok=True)
    npy_path = dccm_dir / f"{safe}_rep{replicate:02d}_dccm.npy"
    np.save(str(npy_path), dccm)
    logging.info(f"  Saved DCCM → {npy_path}")

    # Map domain residue lists to Cα indices (accounting for resid offset)
    def _residues_to_indices(target_positions: list[int]) -> list[int]:
        target_resids = {p + resid_offset for p in target_positions}
        return [i for i, r in enumerate(resids) if r in target_resids]

    idx_nnbp   = _residues_to_indices(_NNBP_RESIDUES)
    idx_fingers = _residues_to_indices(_FINGERS_RESIDUES)
    idx_palm   = _residues_to_indices(_PALM_RESIDUES)
    idx_thumb  = _residues_to_indices(_THUMB_RESIDUES)

    coupling_fingers = _mean_coupling(dccm, idx_nnbp, idx_fingers)
    coupling_palm    = _mean_coupling(dccm, idx_nnbp, idx_palm)
    coupling_thumb   = _mean_coupling(dccm, idx_nnbp, idx_thumb)

    return {
        "mutation": mutation,
        "safe_label": safe,
        "replicate": replicate,
        "n_frames": len(pos_list),
        "n_ca": int(ca.n_atoms),
        "nnbp_fingers_coupling": coupling_fingers,
        "nnbp_palm_coupling": coupling_palm,
        "nnbp_thumb_coupling": coupling_thumb,
        "dccm_path": str(npy_path),
    }


def _plot_dccm(npy_path: Path, label: str, plots_dir: Path) -> None:
    import matplotlib.pyplot as plt

    dccm = np.load(str(npy_path))
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(dccm, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto", origin="lower")
    plt.colorbar(im, ax=ax, label="Cross-correlation")
    ax.set_title(f"DCCM: {label}", fontsize=10, fontweight="bold")
    ax.set_xlabel("Residue index (Cα)")
    ax.set_ylabel("Residue index (Cα)")
    fig.tight_layout()
    out_path = plots_dir / f"{label}_dccm_heatmap.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logging.info(f"  Saved DCCM plot → {out_path}")


def _plot_coupling_summary(summary_df: pd.DataFrame, plots_dir: Path) -> None:
    import matplotlib.pyplot as plt

    if summary_df.empty:
        return

    mutations = sorted(
        summary_df["mutation"].unique(),
        key=lambda m: (m != "WT", "+" not in m, m),
    )

    fig, axes = plt.subplots(1, 3, figsize=(14, 4), sharey=False)
    coupling_cols = [
        ("nnbp_fingers_coupling", "NNBP–Fingers"),
        ("nnbp_palm_coupling",    "NNBP–Palm"),
        ("nnbp_thumb_coupling",   "NNBP–Thumb"),
    ]

    for ax, (col, title) in zip(axes, coupling_cols):
        means = []
        errs = []
        for mut in mutations:
            vals = summary_df[summary_df["mutation"] == mut][col].dropna()
            means.append(float(vals.mean()) if len(vals) else float("nan"))
            errs.append(float(vals.std()) if len(vals) > 1 else 0.0)

        colors = ["#2196F3" if m == "WT" else "#FF5722" for m in mutations]
        x = np.arange(len(mutations))
        bars = ax.bar(x, means, yerr=errs, color=colors, capsize=3, alpha=0.8)

        ax.set_xticks(x)
        ax.set_xticklabels(mutations, rotation=45, ha="right", fontsize=8)
        ax.set_title(title, fontweight="bold")
        ax.set_ylabel("Mean |DCCM| (0–1)")
        ax.grid(axis="y", alpha=0.3, linestyle=":")

    fig.suptitle(
        "Allosteric Coupling of NNBP to RT Polymerase Domains",
        fontsize=12, fontweight="bold",
    )
    fig.tight_layout()
    out_path = plots_dir / "nnbp_allosteric_coupling_summary.png"
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    import matplotlib.pyplot as plt_
    plt_.close(fig)
    logging.info(f"Wrote {out_path}")


def main() -> int:
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    parser = argparse.ArgumentParser(
        description="Compute DCCM and NNBP–domain allosteric coupling scores."
    )
    parser.add_argument("--manifest", type=Path, default=Path("manifests/md_manifest.csv"))
    parser.add_argument("--resid-offset", type=int, default=-3)
    parser.add_argument("--frame-stride", type=int, default=2,
                        help="Frame stride (default 2 = every other frame for speed)")
    parser.add_argument(
        "--mutations",
        nargs="*",
        default=None,
        help="Subset of mutations to process (default: all)",
    )
    parser.add_argument(
        "--dccm-dir",
        type=Path,
        default=Path("results/dccm"),
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=Path("results/dccm_allosteric_coupling.csv"),
    )
    parser.add_argument("--plots-dir", type=Path, default=Path("results/plots/dccm"))
    parser.add_argument("--save-plots", action="store_true",
                        help="Generate per-replicate DCCM heatmap images (slow)")
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

    results: list[dict] = []
    for _, row in mf.iterrows():
        mut = str(row["mutation"])
        rep = int(row["replicate"])
        logging.info(f"Processing {mut} rep{rep} (DCCM)...")
        try:
            result = _process_replicate(
                row, repo_root, args.resid_offset, args.frame_stride, args.dccm_dir
            )
            results.append(result)
            logging.info(
                f"  coupling: fingers={result['nnbp_fingers_coupling']:.3f} "
                f"palm={result['nnbp_palm_coupling']:.3f} "
                f"thumb={result['nnbp_thumb_coupling']:.3f}"
            )
            if args.save_plots:
                safe = str(row.get("safe_label", mut.replace("+", "_")))
                _plot_dccm(
                    Path(result["dccm_path"]),
                    f"{safe}_rep{rep:02d}",
                    args.plots_dir,
                )
        except Exception as exc:
            logging.error(f"  FAILED {mut} rep{rep}: {exc}")

    if not results:
        logging.error("No DCCM results collected.")
        return 1

    out_df = pd.DataFrame(results)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.output_csv, index=False)
    logging.info(f"Wrote {args.output_csv} ({len(out_df)} rows)")

    args.plots_dir.mkdir(parents=True, exist_ok=True)
    _plot_coupling_summary(out_df, args.plots_dir)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
