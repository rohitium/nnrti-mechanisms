#!/usr/bin/env python3
"""Abstract (genotype-free) versions of the pmx NEQ protocol figures.

The figures in ``protocol/<genotype>/`` are data-driven: they carry that
genotype's numbers and are the record of a specific result. These are the
*method* figures — same physics, no numbers — for explaining the protocol
independently of any one leg.

Three panels (the λ-profile figure has no schematic analogue and is skipped):

  01  thermodynamic cycle, WT -> Mut, with the closure relation
  02  hybrid topology at λ = 0 / λ = 1 (drawn from the real V106A hybrid.pdb,
      since a schematic of dual-topology sticks is less clear than the sticks)
  03  Crooks work distributions for holo and apo, with the three estimators
      given as equations rather than fitted values

Usage:
    PYTHONPATH=. python scripts/fep_pmx/plot_protocol_schematic.py
"""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch

HERE = Path(__file__).resolve().parent
OUT_DIR = Path("results/analysis/fep_pmx/protocol_schematic")

# Estimator definitions, shown in the legend of figure 03 instead of fitted values.
# CGI is spelled out: pmx reports it as "Crooks Gaussian Intersection" -- it fits a
# Gaussian to each work distribution and takes dG as their crossing, which follows
# from the Crooks relation being unity at W = dG.
EST_CGI = r"CGI (Crooks Gaussian Intersection): $P(W_f) \cap P(W_r)$"
EST_BAR = r"BAR (Bennett Acceptance Ratio): $P(W_f)/P(W_r) = e^{\beta(W - \Delta G)}$"
EST_JARZ = r"Jarzynski: $\Delta G = -kT\ln\langle e^{-\beta W_f}\rangle$"


def _load_protocol_module():
    """Reuse the stick-drawing helpers from plot_protocol_figures (not a package)."""
    spec = importlib.util.spec_from_file_location(
        "_fep_protocol_figures", HERE / "plot_protocol_figures.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def plot_cycle(out_path: Path) -> None:
    """WT -> Mut cycle with no numbers, closing on the standard relation.

    Cycle closure, read around the square:
        dG_bind^WT + dG_holo = dG_apo + dG_bind^Mut
    hence
        ddG_bind = dG_holo - dG_apo = dG_bind^Mut - dG_bind^WT
    """
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 7)
    ax.axis("off")
    boxes = [
        (1.1, 4.6, "WT RT · DOR\n(holo, λ = 0)"),
        (6.3, 4.6, "Mut RT · DOR\n(holo, λ = 1)"),
        (1.1, 1.15, "WT RT\n(apo, λ = 0)"),
        (6.3, 1.15, "Mut RT\n(apo, λ = 1)"),
    ]
    for x, y, text in boxes:
        ax.add_patch(
            FancyBboxPatch(
                (x, y), 2.6, 1.45, boxstyle="round,pad=0.08,rounding_size=0.15",
                facecolor="#eef3f8", edgecolor="#1f4e79", lw=1.4,
            )
        )
        ax.text(x + 1.3, y + 0.72, text, ha="center", va="center", fontsize=12)

    def arrow(a, b, label):
        ax.annotate("", xy=b, xytext=a,
                    arrowprops=dict(arrowstyle="-|>", color="#1f4e79", lw=1.6))
        ax.text((a[0] + b[0]) / 2, (a[1] + b[1]) / 2 + 0.18, label,
                ha="center", va="bottom", fontsize=12, color="#1f4e79", fontweight="bold")

    arrow((3.7, 5.35), (6.3, 5.35), r"$\Delta G_{\mathrm{holo}}$")
    arrow((3.7, 1.9), (6.3, 1.9), r"$\Delta G_{\mathrm{apo}}$")

    ax.annotate("", xy=(2.4, 4.6), xytext=(2.4, 2.6),
                arrowprops=dict(arrowstyle="-|>", color="0.35", lw=1.2))
    ax.text(1.05, 3.55, r"$\Delta G_{\mathrm{bind}}^{\mathrm{WT}}$",
            fontsize=12, color="0.35", rotation=90, va="center")
    ax.annotate("", xy=(7.6, 4.6), xytext=(7.6, 2.6),
                arrowprops=dict(arrowstyle="-|>", color="0.35", lw=1.2))
    ax.text(8.85, 3.55, r"$\Delta G_{\mathrm{bind}}^{\mathrm{Mut}}$",
            fontsize=12, color="0.35", rotation=90, va="center")

    ax.set_title("Thermodynamic cycle", fontweight="bold", fontsize=15)
    ax.text(
        5.0, 0.3,
        r"$\Delta\Delta G_{\mathrm{bind}} = \Delta G_{\mathrm{holo}} - \Delta G_{\mathrm{apo}}"
        r" = \Delta G_{\mathrm{bind}}^{\mathrm{Mut}} - \Delta G_{\mathrm{bind}}^{\mathrm{WT}}$",
        ha="center", va="center", fontsize=14,
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#fff6e5", edgecolor="#c47b16"),
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_hybrid(out_path: Path, mod) -> None:
    """λ = 0 / λ = 1 dual topology, larger labels, no explanatory subtitles."""
    hybrid_pdb = mod._hybrid_pdb("wt_to_V106A")
    if hybrid_pdb is None:
        raise FileNotFoundError("wt_to_V106A hybrid.pdb not available locally")
    hy = mod._parse_pdb(hybrid_pdb)
    st = mod._parse_pdb(mod.START_PDB)
    mut = [a for a in hy if a["chain"] == "A" and a["resid"] == 103]
    lig = [a for a in st if a["resn"] == "2KW" and a["name"][0] != "H"]
    lookup = {a["name"]: a["xyz"] for a in mut}
    origin = lookup["CB"]
    v1, v2 = lookup["CG1"] - origin, lookup["CG2"] - origin
    xhat = v1 / np.linalg.norm(v1)
    ytmp = v2 - xhat * np.dot(v2, xhat)
    yhat = ytmp / np.linalg.norm(ytmp)
    xy = {a["name"]: mod._project(a["xyz"][None, :], origin, xhat, yhat)[0] for a in mut}
    xyz_lig = np.vstack([a["xyz"] for a in lig])
    p_lig = mod._project(xyz_lig, origin, xhat, yhat)

    backbone = [("N", "CA"), ("CA", "C"), ("C", "O"), ("CA", "CB")]
    val_bonds = [("CB", "CG1"), ("CB", "CG2")]
    ala_bonds = [("CB", "HV1"), ("CB", "HV2")]

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.8), sharex=True, sharey=True)
    panels = [
        {"ax": axes[0], "title": r"$\lambda = 0$",
         "show_val": True, "ghost_val": False, "show_ala": False, "ghost_ala": True},
        {"ax": axes[1], "title": r"$\lambda = 1$",
         "show_val": False, "ghost_val": True, "show_ala": True, "ghost_ala": False},
    ]
    for panel in panels:
        ax = panel["ax"]
        for i in range(len(lig)):
            for j in range(i + 1, len(lig)):
                if np.linalg.norm(xyz_lig[i] - xyz_lig[j]) < 1.85:
                    ax.plot([p_lig[i, 0], p_lig[j, 0]], [p_lig[i, 1], p_lig[j, 1]],
                            color="#9bb6d3", lw=1.0, zorder=1)
        ax.scatter(p_lig[:, 0], p_lig[:, 1], s=12, c="#9bb6d3", zorder=1)
        mod._draw_sticks(ax, xy, backbone, color="#333333", lw=2.2)
        mod._draw_atoms(ax, xy, ["N", "CA", "C", "O", "CB"], color="#4d4d4d", s=70)
        if panel["show_val"] or panel["ghost_val"]:
            alpha = 1.0 if panel["show_val"] else 0.22
            mod._draw_sticks(ax, xy, val_bonds, color="#2ca02c", lw=2.8, alpha=alpha)
            mod._draw_atoms(ax, xy, ["CG1", "CG2"], color="#2ca02c", s=160, alpha=alpha)
        if panel["show_ala"] or panel["ghost_ala"]:
            alpha = 1.0 if panel["show_ala"] else 0.22
            mod._draw_sticks(ax, xy, ala_bonds, color="#d62728", lw=2.0,
                             ls="--" if panel["ghost_ala"] else "-", alpha=alpha)
            mod._draw_atoms(ax, xy, ["HV1", "HV2"], color="#d62728", s=90, alpha=alpha)
        for name, label in (("CB", "Cβ"), ("CG1", "Cγ1"), ("CG2", "Cγ2"), ("CA", "Cα")):
            if name in xy:
                ax.annotate(label, xy[name], textcoords="offset points", xytext=(5, 5),
                            fontsize=13, color="0.2")
        ax.set_aspect("equal")
        keys = ("CA", "CB", "CG1", "CG2", "HV1", "HV2")
        xs = [xy[n][0] for n in keys if n in xy]
        ys = [xy[n][1] for n in keys if n in xy]
        pad = 2.2
        ax.set_xlim(min(xs) - pad, max(xs) + pad)
        ax.set_ylim(min(ys) - pad, max(ys) + pad)
        ax.annotate("Doravirine", p_lig.mean(axis=0), textcoords="offset points",
                    xytext=(6, 8), fontsize=13, color="#3d6fa0", fontstyle="italic")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(panel["title"], fontweight="bold", fontsize=15)
        for spine in ax.spines.values():
            spine.set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _gauss(x, mu, sigma):
    return np.exp(-0.5 * ((x - mu) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))


def plot_work_distributions(out_path: Path) -> None:
    """Schematic Crooks histograms for holo and apo; estimators as equations.

    Replaces the per-genotype 03_neq_work and 05_crooks_overlap figures. Shapes are
    illustrative: holo is drawn with more forward/reverse overlap than apo, which is
    the qualitative pattern in the real data, but no values are implied.
    """
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.1), sharey=True)
    specs = [
        # Titles sit at the OUTER top corner of each panel so the centred legend
        # cannot cover them.
        {"ax": axes[0], "title": "holo", "tx": 0.02, "ha": "left",
         "dg": r"$\Delta G_{\mathrm{holo}}$",
         "mu_f": 1.05, "mu_r": -1.05, "sig": 1.05},
        {"ax": axes[1], "title": "apo", "tx": 0.98, "ha": "right",
         "dg": r"$\Delta G_{\mathrm{apo}}$",
         "mu_f": 1.5, "mu_r": -1.5, "sig": 0.95},
    ]
    x = np.linspace(-4.6, 4.6, 800)
    for sp in specs:
        ax = sp["ax"]
        yf = _gauss(x, sp["mu_f"], sp["sig"])
        yr = _gauss(x, sp["mu_r"], sp["sig"])
        ax.fill_between(x, yr, color="#e8a87c", alpha=0.55, lw=0, label=r"Reverse Work  $P(W_r)$,   $\lambda: 1 \rightarrow 0$")
        ax.fill_between(x, yf, color="#8fb8de", alpha=0.65, lw=0, label=r"Forward Work  $P(W_f)$,   $\lambda: 0 \rightarrow 1$")
        ax.plot(x, yr, color="#d1651a", lw=2.0)
        ax.plot(x, yf, color="#1f6fb2", lw=2.0)
        # The three estimators coincide at the crossing when the work is Gaussian;
        # a single marker keeps the schematic honest rather than implying a spread.
        ax.axvline(0.0, color="0.25", lw=1.6, ls="--")
        ax.plot([0.0], [_gauss(0.0, sp["mu_f"], sp["sig"])], "o", color="#5b2d8e",
                markersize=8, zorder=5)
        ax.annotate(sp["dg"], (0.0, _gauss(0.0, sp["mu_f"], sp["sig"])),
                    textcoords="offset points", xytext=(8, 12), fontsize=14,
                    color="0.2", fontweight="bold")
        # Headroom so the centred legend clears the curves instead of needing a
        # whitespace band under the figure.
        ax.set_ylim(0, _gauss(sp["mu_f"], sp["mu_f"], sp["sig"]) * 1.72)
        ax.text(sp["tx"], 0.97, sp["title"], transform=ax.transAxes,
                ha=sp["ha"], va="top", fontweight="bold", fontsize=15)
        ax.set_xlabel("work $W$", fontsize=14)
        ax.set_xticks([]); ax.set_yticks([])
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[0].set_ylabel("Probability density $P(W)$", fontsize=14)

    handles, labels = axes[0].get_legend_handles_labels()
    blank = plt.Line2D([], [], linestyle="none")
    handles += [blank, blank, blank]
    labels += [EST_CGI, EST_BAR, EST_JARZ]
    fig.legend(handles, labels, loc="upper center", ncol=1, frameon=True,
               fontsize=8.5, bbox_to_anchor=(0.5, 0.995))
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    mod = _load_protocol_module()

    plot_cycle(args.out_dir / "01_thermodynamic_cycle.png")
    plot_hybrid(args.out_dir / "02_hybrid_topology.png", mod)
    plot_work_distributions(args.out_dir / "03_work_distributions.png")
    for name in ("01_thermodynamic_cycle", "02_hybrid_topology", "03_work_distributions"):
        print(f"Wrote {args.out_dir / (name + '.png')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
