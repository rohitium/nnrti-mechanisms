#!/usr/bin/env python3
"""Supplementary Figures 3-5, which existed as captions with no figure behind them.

  Figure 3  convergence diagnostics: block averages and autocorrelation
  Figure 4  protein-level changes: Ca RMSF and secondary structure, WT-referenced
  Figure 5  unliganded and bound NNIBP volume distributions

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.plot_supplementary_figures
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
RES = ROOT / "results/analysis"
#: unliganded crystal structures, for the reference line in Figure 5
APO_CRYSTAL = (101.7, 103.8)
#: Block lengths in frames. Capped so that every point still averages at
#: least MIN_BLOCKS blocks: with 180 frames a block of 45 leaves only four,
#: and the standard error of four numbers is too noisy to read a plateau off.
MIN_BLOCKS = 10
BLOCKS = [1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 18]

METRIC_LABEL = {
    "chl_ring_burial": "chlorocyanophenyl burial",
    "dor_to_S105_mindist": "DOR–Ser105 distance",
    "dor_to_res106_mindist": "DOR–residue 106 distance",
    "res103_bb_to_triazN": "Lys103 backbone–triazolinone N",
    "res179_to_dor_mindist": "Val179–DOR distance",
    "y188_interplanar_deg": "Tyr188 interplanar angle",
}


def _style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 9,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.dpi": 300, "savefig.bbox": "tight",
    })
    return plt


#: per-frame observables worth a convergence diagnostic, with display labels
TRACE_METRICS = {
    "ca_rmsd_angstrom": "C\u03b1 RMSD",
    "dor_rmsd_angstrom": "DOR RMSD",
    "com_distance_angstrom": "DOR\u2013RT COM distance",
    "ser105_min_angstrom": "DOR\u2013Ser105 distance",
}


def autocorrelation(x: np.ndarray) -> np.ndarray:
    """Normalised autocorrelation function, C(0) = 1."""
    x = np.asarray(x, dtype=float)
    x = x - x.mean()
    if np.allclose(x, 0):
        return np.full(x.size, np.nan)
    n = x.size
    # FFT route; the series are short but this keeps it exact and quick
    f = np.fft.rfft(x, 2 * n)
    acf = np.fft.irfft(f * np.conjugate(f), 2 * n)[:n].real
    return acf / acf[0]


def block_sem(x: np.ndarray, block: int) -> float:
    """Standard error from the means of non-overlapping blocks of `block`.

    Flyvbjerg and Petersen: for correlated data this grows with block size and
    plateaus once a block is longer than the correlation time. The plateau is
    the honest standard error; the naive value at block = 1 is an
    underestimate by the square root of the statistical inefficiency.
    """
    x = np.asarray(x, dtype=float)
    m = x.size // block
    if m < 2:
        return np.nan
    means = x[: m * block].reshape(m, block).mean(axis=1)
    return float(means.std(ddof=1) / np.sqrt(m))


def figure3(plt, out: Path) -> None:
    tr = pd.read_csv(RES / "md_convergence/tables/frame_traces.csv")
    dt = float(tr.total_ns.median() / tr.n_frames.median())      # ns per frame

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.0))
    colours = ["#4C72B0", "#DD8452", "#55A868", "#C44E52"]

    n_series = 0
    for (col, label), c in zip(TRACE_METRICS.items(), colours):
        acfs, ratios = [], []
        for _, s in tr.groupby(["safe_label", "replicate"]):
            x = s.sort_values("frame_index")[col].to_numpy(dtype=float)
            if x.size < 40 or np.allclose(x, x[0]):
                continue
            acfs.append(autocorrelation(x))
            base = block_sem(x, 1)
            usable = [b for b in BLOCKS if x.size // b >= MIN_BLOCKS]
            if base and np.isfinite(base) and len(usable) == len(BLOCKS):
                ratios.append([block_sem(x, b) / base for b in BLOCKS])
        if not acfs:
            continue
        n_series = max(n_series, len(acfs))
        L = min(len(a) for a in acfs)
        A = np.vstack([a[:L] for a in acfs])
        lag = np.arange(L) * dt
        keep = lag <= 30
        ax1.plot(lag[keep], np.nanmean(A, axis=0)[keep], color=c, lw=1.3, label=label)
        ax1.fill_between(lag[keep],
                         np.nanpercentile(A, 25, axis=0)[keep],
                         np.nanpercentile(A, 75, axis=0)[keep],
                         color=c, alpha=0.15, lw=0)

        R = np.vstack([r for r in ratios if np.all(np.isfinite(r))])
        ax2.plot(np.array(BLOCKS) * dt, R.mean(axis=0), color=c, lw=1.3,
                 marker="o", ms=2.5, label=label)

    ax1.axhline(0, color="0.6", lw=0.6)
    ax1.set_xlabel("lag (ns)")
    ax1.set_ylabel("autocorrelation C(\u03c4)")
    ax1.set_title(f"A  Autocorrelation, {n_series} replicates", loc="left")
    ax1.legend(frameon=False, fontsize=6.5)

    ax2.axhline(1.0, color="0.6", lw=0.6, ls=":")
    ax2.set_xlabel("block length (ns)")
    ax2.set_ylabel("blocked SEM / naive SEM")
    ax2.set_title("B  Block averaging", loc="left")

    fig.tight_layout()
    fig.savefig(out / "supplementary_figure_3_convergence.png")
    plt.close(fig)
    print(f"  Figure 3: {n_series} replicates, {len(TRACE_METRICS)} observables")


def figure4(plt, out: Path) -> None:
    rmsf = pd.read_csv(RES / "protein_level/tables/ca_rmsf_per_residue.csv")
    ss = pd.read_csv(RES / "protein_level/tables/secondary_structure_per_residue.csv")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.2, 4.2), sharex=True)

    def wt_ref(d, col):
        m = d.groupby(["safe_label", "canonical"])[col].mean().unstack(0)
        if "wt" not in m.columns:
            raise SystemExit(f"no wild type column in {col}")
        return m.drop(columns="wt").sub(m["wt"], axis=0)

    dr = wt_ref(rmsf, "rmsf")
    ax1.axhline(0, color="0.6", lw=0.6)
    for c in dr.columns:
        ax1.plot(dr.index, dr[c], lw=0.5, alpha=0.4, color="#4C72B0")
    ax1.plot(dr.index, dr.mean(axis=1), lw=1.3, color="#C44E52",
             label="panel mean")
    ax1.set_ylabel("Δ Cα RMSF (Å)")
    ax1.set_title("A  Per-residue flexibility, wild-type referenced", loc="left")
    ax1.legend(frameon=False, loc="upper right")

    ss["ordered"] = ss.frac_helix + ss.frac_strand
    ds = wt_ref(ss, "ordered")
    ax2.axhline(0, color="0.6", lw=0.6)
    for c in ds.columns:
        ax2.plot(ds.index, ds[c], lw=0.5, alpha=0.4, color="#55A868")
    ax2.plot(ds.index, ds.mean(axis=1), lw=1.3, color="#C44E52")
    ax2.set_xlabel("residue (canonical HIV-1 RT numbering)")
    ax2.set_ylabel("Δ ordered fraction")
    ax2.set_title("B  Secondary structure, wild-type referenced", loc="left")

    fig.tight_layout()
    fig.savefig(out / "supplementary_figure_4_protein_level.png")
    plt.close(fig)
    print(f"  Figure 4: {dr.shape[1]} genotypes, {dr.shape[0]} residues")


def figure5(plt, out: Path) -> None:
    holo = pd.read_csv(ROOT / "results/.checkpoints/.checkpoint_pocket_volume_profiles.csv")
    apo = pd.read_csv(RES / "apo_ensembles/tables/apo_frame_metrics.csv")
    apo = apo[apo.window == "matched_10ns"]

    # A handful of frames return the grid-saturation value, where the measuring
    # grid has landed in bulk solvent instead of the pocket. The Calpha cage
    # cannot enclose more than about 300 A^3, so anything above MAX_PHYSICAL is
    # not a measurement. F227C replicate 2 supplies 9 of these, reaching
    # 4172 A^3, which on its own sets the y axis and hides every distribution.
    MAX_PHYSICAL = 500.0
    n_h = int((holo.pocket_volume_proxy_angstrom3 > MAX_PHYSICAL).sum())
    n_a = int((apo.pocket_volume > MAX_PHYSICAL).sum())
    if n_h or n_a:
        print(f"  Figure 5: excluded {n_h} holo and {n_a} apo frames above "
              f"{MAX_PHYSICAL:.0f} A^3 as grid-saturation artefacts")
    holo = holo[holo.pocket_volume_proxy_angstrom3 <= MAX_PHYSICAL]
    apo = apo[apo.pocket_volume <= MAX_PHYSICAL]

    labels = sorted(apo.safe_label.unique())
    fig, ax = plt.subplots(figsize=(7.2, 3.2))
    pos = np.arange(len(labels))
    hol, ap = [], []
    for lab in labels:
        h = holo[holo.safe_label.str.lower() == lab.lower()].pocket_volume_proxy_angstrom3
        hol.append(h.to_numpy() if len(h) else np.array([np.nan]))
        ap.append(apo[apo.safe_label == lab].pocket_volume.to_numpy())

    for data, off, color, name in ((hol, -0.19, "#4C72B0", "bound (holo)"),
                                   (ap, 0.19, "#DD8452", "unliganded (apo)")):
        parts = ax.violinplot([d[~np.isnan(d)] for d in data],
                              positions=pos + off, widths=0.34,
                              showmedians=True, showextrema=False)
        for b in parts["bodies"]:
            b.set_facecolor(color); b.set_alpha(0.75); b.set_edgecolor("none")
        parts["cmedians"].set_color("0.2"); parts["cmedians"].set_linewidth(0.9)
        ax.plot([], [], color=color, lw=5, alpha=0.75, label=name)

    ax.axhspan(APO_CRYSTAL[0], APO_CRYSTAL[1], color="0.5", alpha=0.25, zorder=0)
    ax.text(len(labels) - 0.5, APO_CRYSTAL[1] + 4, "unliganded crystal (1DLO, 1HMV)",
            ha="right", fontsize=6.5, color="0.35")
    ax.set_xticks(pos)
    ax.set_xticklabels([l.replace("_", "+").upper() if l != "wt" else "WT"
                        for l in labels], rotation=30, ha="right")
    ax.set_ylabel("NNIBP volume (Å³)")
    ax.set_title("Unliganded and bound NNIBP volume", loc="left")
    ax.legend(frameon=False, loc="upper right", ncol=2)
    fig.tight_layout()
    fig.savefig(out / "supplementary_figure_5_pocket_volume.png")
    plt.close(fig)
    print(f"  Figure 5: {len(labels)} genotypes with apo data")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--output-dir", type=Path,
                    default=RES / "supplementary_figures")
    args = ap.parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plt = _style()
    figure3(plt, args.output_dir)
    figure4(plt, args.output_dir)
    figure5(plt, args.output_dir)
    print(f"wrote to {args.output_dir.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
