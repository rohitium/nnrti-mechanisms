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
#: shared lag axis in ns for the autocorrelation panels
LAG_GRID = np.arange(0.0, 30.01, 0.25)
#: shared time axis in ns for the raw traces
TIME_GRID = np.arange(0.0, 100.01, 0.5)
#: Block lengths in frames. Block means are pooled across a genotype's three
#: replicates, so a block of n frames yields 3 * (N // n) blocks rather than
#: N // n, which is what makes the longer blocks usable at all: the
#: autocorrelation of Ca RMSD runs to 20-30 ns, and a plateau can only appear
#: once a block is longer than that. At 60 frames (~25 ns) this still leaves
#: only 9 blocks, so the far end of the curve stays noisy by construction --
#: the trajectories are not long enough for it to be otherwise.
#: Block lengths in NANOSECONDS, not frames. The series do not share a frame
#: rate -- the same 100 ns appears as 139, 180, 238, 360 or 540 frames -- so a
#: fixed frame count would mean 50 ns for one replicate and 25 ns for another,
#: and the curve would mix durations. Each series converts these to its own
#: frame count. At 50 ns a 100 ns replicate yields 2 blocks, so 6 across the
#: three replicates, which is the longest block the data can support.
MIN_BLOCKS = 6
BLOCKS_NS = [0.5, 1, 2, 3, 5, 7, 10, 14, 20, 25, 33, 40, 50]

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


#: the three observables shown per-genotype: ligand-internal, ligand-to-protein
#: and protein-global. Ser105 distance behaves between the first two.
PANEL_METRICS = {
    "dor_rmsd_angstrom": "DOR RMSD",
    "com_distance_angstrom": "DOR\u2013RT COM distance",
    "ca_rmsd_angstrom": "C\u03b1 RMSD",
}


def figure1(plt, out: Path) -> None:
    """Nine panels: traces, autocorrelation, statistical inefficiency."""
    tr = pd.read_csv(RES / "md_convergence/tables/frame_traces.csv")

    fig, axes = plt.subplots(3, 3, figsize=(7.2, 6.4))
    for j, (col, label) in enumerate(PANEL_METRICS.items()):
        a_raw, a_acf, a_g = axes[0, j], axes[1, j], axes[2, j]
        acf_by_g, g_by_g, raw_by_g, block_counts = {}, {}, {}, []

        for gen, sub in tr.groupby("safe_label"):
            acfs, series, dts, traces = [], [], [], []
            for _, r in sub.groupby("replicate"):
                r = r.sort_values("frame_index")
                x = r[col].to_numpy(dtype=float)
                if x.size < 40 or np.allclose(x, x[0]):
                    continue
                d = float(r.total_ns.iloc[0]) / float(x.size)
                a = autocorrelation(x)
                acfs.append(np.interp(LAG_GRID, np.arange(a.size) * d, a,
                                      left=np.nan, right=np.nan))
                traces.append(np.interp(TIME_GRID, np.arange(x.size) * d, x,
                                        left=np.nan, right=np.nan))
                series.append(x)
                dts.append(d)
            if not series:
                continue
            raw_by_g[gen] = np.nanmean(np.vstack(traces), axis=0)
            acf_by_g[gen] = np.nanmean(np.vstack(acfs), axis=0)

            base = np.mean([block_sem(x, 1) for x in series])
            ratio, counts = [], []
            for b_ns in BLOCKS_NS:
                chunks = []
                for x, d in zip(series, dts):
                    b = max(1, int(round(b_ns / d)))
                    if x.size // b < 2:
                        continue
                    chunks.append(
                        x[: (x.size // b) * b].reshape(-1, b).mean(axis=1))
                means = np.concatenate(chunks) if chunks else np.array([])
                counts.append(means.size)
                # statistical inefficiency g = (blocked SEM / naive SEM)^2
                ratio.append((means.std(ddof=1) / np.sqrt(means.size) / base) ** 2
                             if means.size >= MIN_BLOCKS and base else np.nan)
            g_by_g[gen] = np.array(ratio)
            block_counts.append(counts)

        R = np.vstack(list(raw_by_g.values()))
        for row in R:
            a_raw.plot(TIME_GRID, row, color="#937860", lw=0.5, alpha=0.35)
        a_raw.plot(TIME_GRID, np.nanmedian(R, axis=0), color="#C44E52", lw=1.4)
        a_raw.set_title(label, loc="left", fontsize=8.5)
        a_raw.set_xlabel("time (ns)")
        if j == 0:
            a_raw.set_ylabel("value (\u00c5)")

        M = np.vstack(list(acf_by_g.values()))
        for row in M:
            a_acf.plot(LAG_GRID, row, color="#4C72B0", lw=0.5, alpha=0.35)
        a_acf.plot(LAG_GRID, np.nanmedian(M, axis=0), color="#C44E52", lw=1.4)
        a_acf.axhline(0, color="0.6", lw=0.6)
        a_acf.set_xlabel("lag \u03c4 (ns)")
        if j == 0:
            a_acf.set_ylabel("C(\u03c4)")

        bl = np.array(BLOCKS_NS)
        n_blocks = np.nanmedian(np.vstack(block_counts), axis=0)
        B = np.vstack(list(g_by_g.values()))
        for row in B:
            a_g.plot(bl, row, color="#55A868", lw=0.5, alpha=0.3)
        med = np.nanmedian(B, axis=0)
        a_g.plot(bl, med, color="#C44E52", lw=1.4)
        from scipy import stats as _st
        m_arr = np.asarray(n_blocks, dtype=float)
        # g is the square of a ratio of standard deviations, so the chi-square
        # limits on the standard deviation are squared too
        lo = (m_arr - 1) / _st.chi2.ppf(0.975, m_arr - 1)
        hi = (m_arr - 1) / _st.chi2.ppf(0.025, m_arr - 1)
        a_g.fill_between(bl, med * lo, med * hi, color="#C44E52",
                         alpha=0.16, lw=0)
        a_g.axhline(1.0, color="0.6", lw=0.6, ls=":")
        a_g.set_xlabel("block length (ns)")
        if j == 0:
            a_g.set_ylabel("statistical inefficiency $g$")
        a_g.text(bl[-1], a_g.get_ylim()[1] * 0.96, f"{int(n_blocks[-1])} blocks",
                 fontsize=5.5, color="0.4", va="top", ha="right")

    axes[0, 0].plot([], [], color="0.5", lw=0.8, alpha=0.6, label="genotype")
    axes[0, 0].plot([], [], color="#C44E52", lw=1.4, label="median")
    axes[0, 0].legend(frameon=False, fontsize=6.5, loc="best")
    for ax, tag in zip(axes.ravel(), "ABCDEFGHI"):
        ax.text(-0.17, 1.04, tag, transform=ax.transAxes,
                fontsize=9, fontweight="bold", va="bottom")
    fig.tight_layout()
    fig.savefig(out / "supplementary_figure_1_convergence.png")
    plt.close(fig)
    print(f"  Figure 1: 9 panels, {len(raw_by_g)} genotypes")


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
    figure1(plt, args.output_dir)
    figure4(plt, args.output_dir)
    figure5(plt, args.output_dir)
    print(f"wrote to {args.output_dir.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
