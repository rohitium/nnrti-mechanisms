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
#: Ca RMSD is deliberately absent. The reviewer asked about the sampled DOR
#: pose, and a protein-global coordinate answers a different question while
#: making the figure harder to read.
PANEL_METRICS = {
    "dor_rmsd_angstrom": "DOR RMSD",
    "com_distance_angstrom": "DOR\u2013RT COM distance",
}


def _fold_map(xlsx: Path) -> dict:
    t = pd.read_excel(xlsx, header=1)
    t.columns = ["mutation", "rpv_fold", "dor_fold"]
    t = t.dropna(subset=["mutation", "dor_fold"])
    k = t["mutation"].astype(str).str.replace(r",\s*", "+", regex=True).str.strip()
    return dict(zip(k, t["dor_fold"].astype(float)))


#: The submitted manuscript's convergence panels, reused as-is. These are read
#: and never rewritten: compositing the published PNGs is the only way to
#: guarantee the top row is identical to what was submitted, rather than a
#: re-plot that happens to look the same.
SUBMITTED_PANELS = [
    ("dor_rmsd_convergence.png", "DOR RMSD"),
    ("com_distance_convergence.png", "DOR\u2013RT COM distance"),
]


def _fold_map(xlsx: Path) -> dict:
    t = pd.read_excel(xlsx, header=1)
    t.columns = ["mutation", "rpv_fold", "dor_fold"]
    t = t.dropna(subset=["mutation", "dor_fold"])
    k = t["mutation"].astype(str).str.replace(r",\s*", "+", regex=True).str.strip()
    return dict(zip(k, t["dor_fold"].astype(float)))


#: Exactly the settings plot_convergence_panel used for the submitted panels.
#: The new panels are rendered as their own 10 x 6 figures with these values
#: and then composited, so every panel in the finished figure carries the same
#: typeface at the same physical size. Matching by eye in a shared grid cannot
#: do that, because the submitted panels are images and scale with the cell.
PANEL_FIGSIZE = (10.0, 6.0)
PANEL_LABEL_FS = 14
PANEL_TICK_FS = 12
PANEL_LEGEND_FS = 11
PANEL_TITLE_FS = 15


def _fold_map(xlsx: Path) -> dict:
    t = pd.read_excel(xlsx, header=1)
    t.columns = ["mutation", "rpv_fold", "dor_fold"]
    t = t.dropna(subset=["mutation", "dor_fold"])
    k = t["mutation"].astype(str).str.replace(r",\s*", "+", regex=True).str.strip()
    return dict(zip(k, t["dor_fold"].astype(float)))


def _series_stats(tr, col):
    """Per-genotype autocorrelation and blocked standard error."""
    out = {}
    for gen, sub in tr.groupby("mutation"):
        acfs, series, dts = [], [], []
        for _, r in sub.groupby("replicate"):
            r = r.sort_values("frame_index")
            x = r[col].to_numpy(dtype=float)
            if x.size < 40 or np.allclose(x, x[0]):
                continue
            d = float(r.total_ns.iloc[0]) / float(x.size)
            a = autocorrelation(x)
            acfs.append(np.interp(LAG_GRID, np.arange(a.size) * d, a,
                                  left=np.nan, right=np.nan))
            series.append(x)
            dts.append(d)
        if not series:
            continue
        A = np.vstack(acfs)
        se = []
        for b_ns in BLOCKS_NS:
            chunks = []
            for x, d in zip(series, dts):
                b = max(1, int(round(b_ns / d)))
                if x.size // b < 2:
                    continue
                chunks.append(x[: (x.size // b) * b].reshape(-1, b).mean(axis=1))
            m = np.concatenate(chunks) if chunks else np.array([])
            se.append(m.std(ddof=1) / np.sqrt(m.size)
                      if m.size >= MIN_BLOCKS else np.nan)
        out[gen] = (np.nanmean(A, 0),
                    np.nanstd(A, 0, ddof=1) / np.sqrt(A.shape[0]), np.array(se))
    return out


def _render_panel(plt, path, kind, stats, mutants, colour, title, ylabel):
    WT_C = "#B2182B"
    fig, ax = plt.subplots(figsize=PANEL_FIGSIZE)
    x = LAG_GRID if kind == "acf" else np.array(BLOCKS_NS)
    for gen in mutants:
        if gen not in stats:
            continue
        mean, sem, se = stats[gen]
        c = colour(gen)
        if kind == "acf":
            ax.fill_between(x, mean - sem, mean + sem, color=c, alpha=0.10,
                            lw=0, zorder=2)
            ax.plot(x, mean, color=c, lw=1.4, zorder=3)
        else:
            ax.plot(x, se, color=c, lw=1.4, marker="o", ms=4, zorder=3)
    if "WT" in stats:
        mean, sem, se = stats["WT"]
        if kind == "acf":
            ax.fill_between(x, mean - sem, mean + sem, color=WT_C, alpha=0.18,
                            lw=0, zorder=4)
            ax.plot(x, mean, color=WT_C, lw=2.6, zorder=5, label="WT")
        else:
            ax.plot(x, se, color=WT_C, lw=2.6, marker="o", ms=5, zorder=5,
                    label="WT")
    if kind == "acf":
        ax.axhline(0, color="0.5", lw=0.8)
        ax.set_xlim(0, LAG_GRID[-1])
        ax.set_xlabel("Time lag (ns)", fontsize=PANEL_LABEL_FS)
    else:
        ax.set_xlabel("Block length (ns)", fontsize=PANEL_LABEL_FS)
    ax.set_ylabel(ylabel, fontsize=PANEL_LABEL_FS)
    ax.set_title(title, fontsize=PANEL_TITLE_FS, pad=12)
    ax.tick_params(labelsize=PANEL_TICK_FS)
    ax.grid(alpha=0.25, lw=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.legend(fontsize=PANEL_LEGEND_FS, frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def figure1(plt, out: Path) -> None:
    """Submitted traces on top; new panels rendered to match, then composited."""
    from matplotlib.colors import LogNorm
    import matplotlib.image as mpimg

    plot_dir = RES / "md_convergence/plots"
    scratch = out / "_panels"
    scratch.mkdir(parents=True, exist_ok=True)

    tr = pd.read_csv(RES / "md_convergence/tables/frame_traces.csv")
    tr = tr[tr.mutation != "F227C"]
    folds = _fold_map(ROOT / "data/DRM-susceptibilities.csv.xlsx")
    mutants = sorted((m for m in tr.mutation.unique() if m != "WT"),
                     key=lambda m: folds.get(m, np.nan))
    vals = np.array([folds[m] for m in mutants if m in folds], dtype=float)
    norm = LogNorm(vmin=vals.min(), vmax=vals.max())
    ramp = plt.get_cmap("Blues")
    colour = lambda m: (ramp(0.35 + 0.65 * norm(folds[m])) if m in folds else "0.6")

    TITLES = {"dor_rmsd_angstrom": "Autocorrelation of doravirine internal RMSD",
              "com_distance_angstrom": "Autocorrelation of the doravirine\u2013RT separation"}
    SE_TITLES = {"dor_rmsd_angstrom": "Blocked standard error of doravirine internal RMSD",
                 "com_distance_angstrom": "Blocked standard error of the doravirine\u2013RT separation"}
    YL = {"dor_rmsd_angstrom": "DOR RMSD (\u00c5)",
          "com_distance_angstrom": "DOR\u2013RT COM distance (\u00c5)"}

    paths = [plot_dir / SUBMITTED_PANELS[0][0], plot_dir / SUBMITTED_PANELS[1][0]]
    for kind, titles in (("acf", TITLES), ("se", SE_TITLES)):
        for col in PANEL_METRICS:
            stats = _series_stats(tr, col)
            q = scratch / f"{kind}_{col}.png"
            _render_panel(plt, q, kind, stats, mutants, colour, titles[col],
                          "Autocorrelation" if kind == "acf" else YL[col])
            paths.append(q)

    imgs = [mpimg.imread(str(q)) for q in paths]
    ar = imgs[0].shape[0] / imgs[0].shape[1]
    W = 6.0
    fig, axes = plt.subplots(3, 2, figsize=(2 * W, 3 * W * ar))
    for ax, im, tag in zip(axes.ravel(), imgs, "ABCDEF"):
        ax.imshow(im)
        ax.set_axis_off()
        ax.text(-0.02, 1.0, tag, transform=ax.transAxes, fontsize=18,
                fontweight="bold", va="top", ha="right")
    fig.subplots_adjust(left=0.03, right=0.995, top=0.99, bottom=0.005,
                        wspace=0.02, hspace=0.04)
    fig.savefig(out / "supplementary_figure_1_convergence.png", dpi=300)
    plt.close(fig)
    print(f"  Figure 1: 6 panels composited, all rendered at "
          f"{PANEL_FIGSIZE[0]:.0f}x{PANEL_FIGSIZE[1]:.0f} in with matching fonts")


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
