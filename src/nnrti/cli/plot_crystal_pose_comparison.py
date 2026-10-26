#!/usr/bin/env python3
"""Four views of the 4NCG / 7Z2G agreement.

Reviewer 1 doubts that one DOR-bound crystal structure stands for the
biologically relevant RT-DOR ensemble. compare_dor_crystal_poses answers the
narrow version of that doubt -- is the 4NCG pose specific to 4NCG -- and this
plots the answer:

    a  where the two depositions differ, by subdomain, in the palm frame. The
       differences sit in the subdomains that grip nucleic acid, which only one
       of the two structures has. The NNIBP lining is the flattest bar.
    b  ligand-to-residue minimum heavy-atom distances, one point per NNIBP
       residue. Superposition-free, so it cannot be flattered by the fit.
    c  per-residue side-chain RMSD in the NNIBP frame, with the 3.65 A
       coordinate uncertainty of the cryo-EM structure marked.
    d  the four DOR dihedrals: the WT simulated distribution, with both
       experimental values on it.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.plot_crystal_pose_comparison
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REF_C = "#2c7fb8"   # 4NCG, X-ray
MOB_C = "#d7301f"   # 7Z2G, cryo-EM
NNIBP_C = "#238b45"


def main(argv: list[str] | None = None) -> int:
    root = Path(__file__).resolve().parents[3]
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--tables", type=Path,
                   default=root / "results/analysis/crystal_pose_comparison/tables")
    p.add_argument("--torsion-frames", type=Path,
                   default=root / "results/analysis/dor_torsions/tables/"
                                  "dor_torsions_per_frame.csv")
    p.add_argument("--output-dir", type=Path,
                   default=root / "results/analysis/crystal_pose_comparison/plots")
    p.add_argument("--reference", default="4NCG")
    p.add_argument("--mobile", default="7Z2G")
    p.add_argument("--dpi", type=int, default=300)
    args = p.parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    sub = pd.read_csv(args.tables / "subdomain_rmsd.csv")
    res = pd.read_csv(args.tables / "nnibp_residue_comparison.csv")
    tor = pd.read_csv(args.tables / "dor_torsions_crystal.csv")
    lig = pd.read_csv(args.tables / "ligand_rmsd.csv").set_index("frame")

    ref_col = f"min_dist_to_DOR_{args.reference}_A"
    mob_col = f"min_dist_to_DOR_{args.mobile}_A"
    labels = [f"{r.resname.capitalize()}{r.resid}" for r in res.itertuples()]

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.0))
    (ax_a, ax_b), (ax_c, ax_d) = axes

    # -- a: subdomain Ca RMSD ------------------------------------------------ #
    s = sub.copy()
    s["short"] = (s["subdomain"].str.replace("p66 ", "", regex=False)
                  .str.replace(" lining (15 residues)", "\nlining", regex=False))
    is_nnibp = s["subdomain"].str.contains("NNIBP")
    ax_a.bar(s["short"], s["ca_rmsd_palm_frame_A"],
             color=[NNIBP_C if n else "#bdbdbd" for n in is_nnibp],
             edgecolor="black", linewidth=0.5)
    for x, (v, n) in enumerate(zip(s["ca_rmsd_palm_frame_A"], s["n_residues"])):
        ax_a.text(x, v + 0.04, f"{v:.2f}", ha="center", fontsize=7.5)
    ax_a.set_ylabel(r"C$\alpha$ RMSD, palm frame ($\rm\AA$)")
    ax_a.set_title("a   where the two structures differ", loc="left",
                   fontweight="bold", fontsize=10)
    ax_a.tick_params(axis="x", labelsize=8, rotation=30)
    for t in ax_a.get_xticklabels():
        t.set_ha("right")
    ax_a.set_ylim(0, s["ca_rmsd_palm_frame_A"].max() * 1.2)
    ax_a.spines[["top", "right"]].set_visible(False)

    # -- b: contact distances, superposition-free ---------------------------- #
    # a dumbbell rather than a scatter: 16 residues cluster between 3 and 4 A,
    # where point labels would sit on top of one another
    order = res.sort_values(ref_col).reset_index(drop=True)
    y = np.arange(len(order))
    ax_b.hlines(y, order[ref_col], order[mob_col], color="#999999", lw=1.4,
                zorder=1)
    ax_b.scatter(order[ref_col], y, s=38, color=REF_C, edgecolor="black",
                 linewidth=0.5, zorder=3, label=f"{args.reference} (X-ray)")
    ax_b.scatter(order[mob_col], y, s=38, color=MOB_C, edgecolor="black",
                 linewidth=0.5, zorder=3, label=f"{args.mobile} (cryo-EM)")
    ax_b.set_yticks(y)
    ax_b.set_yticklabels([f"{t.resname.capitalize()}{t.resid}"
                          for t in order.itertuples()], fontsize=8)
    ax_b.invert_yaxis()
    r = res[ref_col].corr(res[mob_col])
    mad = (res[mob_col] - res[ref_col]).abs().median()
    biggest = (res[mob_col] - res[ref_col]).abs().max()
    ax_b.set_xlabel(r"minimum heavy-atom distance to DOR ($\rm\AA$)")
    ax_b.set_title("b   closest approach of each residue to DOR", loc="left",
                   fontweight="bold", fontsize=10)
    ax_b.text(0.97, 0.78,
              f"$r$ = {r:.3f}\nmedian |$\\Delta$| = {mad:.2f} $\\rm\\AA$"
              f"\nlargest |$\\Delta$| = {biggest:.2f} $\\rm\\AA$",
              transform=ax_b.transAxes, va="top", ha="right", fontsize=8.5)
    ax_b.legend(loc="upper right", bbox_to_anchor=(1.0, 1.02), fontsize=8,
                frameon=False)
    ax_b.spines[["top", "right"]].set_visible(False)

    # -- c: per-residue side-chain RMSD -------------------------------------- #
    sc = res.dropna(subset=["sidechain_rmsd_A"])
    sc_labels = [f"{t.resname.capitalize()}{t.resid}" for t in sc.itertuples()]
    ax_c.bar(sc_labels, sc["sidechain_rmsd_A"], color="#9ecae1",
             edgecolor="black", linewidth=0.5)
    ax_c.axhline(sc["sidechain_rmsd_A"].median(), color=MOB_C, lw=1.2, ls="--",
                 label=f"median {sc['sidechain_rmsd_A'].median():.2f} $\\rm\\AA$")
    ax_c.set_ylabel(r"side-chain RMSD, NNIBP frame ($\rm\AA$)")
    ax_c.set_title("c   NNIBP side chains", loc="left", fontweight="bold",
                   fontsize=10)
    ax_c.tick_params(axis="x", labelsize=7.5, rotation=60)
    for t in ax_c.get_xticklabels():
        t.set_ha("right")
    ax_c.legend(fontsize=8, frameon=False)
    ax_c.spines[["top", "right"]].set_visible(False)

    # -- d: DOR dihedrals against the WT ensemble ---------------------------- #
    per_frame = pd.read_csv(args.torsion_frames)
    wt = per_frame[per_frame["safe_label"].str.lower() == "wt"]
    names = list(tor["torsion"])
    parts = ax_d.violinplot([wt[n].to_numpy() for n in names],
                            positions=range(len(names)), widths=0.8,
                            showextrema=False, showmedians=False)
    for body in parts["bodies"]:
        body.set_facecolor("#cccccc")
        body.set_edgecolor("#888888")
        body.set_alpha(0.9)
    for i, row in tor.iterrows():
        ax_d.plot([i - 0.33, i + 0.33], [row[f"{args.reference}_deg"]] * 2,
                  color=REF_C, lw=2.2,
                  label=args.reference if i == 0 else None)
        ax_d.plot([i - 0.33, i + 0.33], [row[f"{args.mobile}_deg"]] * 2,
                  color=MOB_C, lw=2.2, ls="--",
                  label=args.mobile if i == 0 else None)
    ax_d.set_xticks(range(len(names)))
    ax_d.set_xticklabels([rf"$\tau_{n[-1]}$" for n in names])
    ax_d.set_ylabel("dihedral (deg)")
    ax_d.set_title(f"d   DOR dihedrals: WT simulated ensemble "
                   f"(n = {len(wt)} frames)", loc="left", fontweight="bold",
                   fontsize=10)
    ax_d.legend(fontsize=8, frameon=False, loc="lower right")
    ax_d.spines[["top", "right"]].set_visible(False)

    nn = lig.loc["nnibp", "ligand_rmsd_A"]
    fig.suptitle(
        f"{args.reference} (X-ray, 2.58 $\\rm\\AA$) vs {args.mobile} "
        f"(cryo-EM, 3.65 $\\rm\\AA$, + DNA): DOR heavy-atom RMSD "
        f"{nn:.2f} $\\rm\\AA$ in the NNIBP frame",
        fontsize=11, y=0.985)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    for ext in ("png", "pdf"):
        fig.savefig(args.output_dir / f"crystal_pose_comparison.{ext}",
                    dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"written: {args.output_dir}/crystal_pose_comparison.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
