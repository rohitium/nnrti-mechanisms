#!/usr/bin/env python3
"""Replot / reinterpret modern_md_suite outputs (no traj reload).

Fixes:
  - ligand RMSF from saved npy (md.rmsf was wrong — C1x artifact)
  - H-bond Δ heatmap excluding self-mutation sites; K103 confirmation bars
  - pocket volume vs experimental fold scatter
  - frame-level DOR conformational map (PCA density per genotype)
  - clearer contact-network view for story mutants
  - DCCM one-pager with plain-language caption

    python -m src.analysis.cli.replot_modern_md_suite
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scripts.fep_pmx.combine_neq import load_experimental, EXPERIMENTAL_CSV

REPO = Path(__file__).resolve().parents[3]
OUT = Path("results/analysis/modern_md_suite")
NNIBP_AUTH = (100, 101, 103, 106, 107, 108, 179, 181, 188, 189, 190, 227, 229, 234, 318)
STORY = ("WT", "V106A", "V106I", "V106I+F227C", "V106A+F227L", "F227C", "A98G+F227C", "G190E", "Y188L")

# Genotypes that mutate a given auth residue — Δoccupancy at that site is mostly chemical change.
MUTATED_SITE = {
    "F227C": {227},
    "A98G+F227C": {98, 227},
    "V106I+F227C": {106, 227},
    "V106A": {106},
    "V106I": {106},
    "V106M": {106},
    "V106A+F227L": {106, 227},
    "V106A+L234I": {106, 234},
    "V106A+P225H": {106, 225},
    "G190A": {190},
    "G190S": {190},
    "G190E": {190},
    "Y181C": {181},
    "Y188L": {188},
    "Y318F": {318},
    "K103N": {103},
    "K103N+M230L": {103, 230},
    "K103N+P225H": {103, 225},
    "L100I+K103N": {100, 103},
}


def _mutation_sort_key(mutation: str) -> tuple[int, str]:
    if mutation == "WT":
        return (0, mutation)
    if "+" in mutation:
        return (2, mutation)
    return (1, mutation)


def recompute_rmsf_from_npy(npy: Path, inv: pd.DataFrame) -> pd.DataFrame:
    """Manual RMSF of pocket-aligned ligand heavy atoms saved during the suite run."""
    name_path = OUT / "tables" / "ligand_rmsf_per_atom.csv"
    names = None
    if name_path.is_file():
        names = (
            pd.read_csv(name_path)
            .query("mutation == 'WT' and replicate == 1")
            .sort_values("atom_index")["atom_name"]
            .tolist()
        )
    rows = []
    for _, r in inv.iterrows():
        lig = np.load(npy / f"{r['npy_stem']}_lig_heavy_xyz.npy")
        mean = lig.mean(axis=0)
        rmsf = np.sqrt(((lig - mean) ** 2).sum(axis=-1).mean(axis=0))
        atom_names = names if names and len(names) == lig.shape[1] else [f"a{i}" for i in range(lig.shape[1])]
        for i, (name, val) in enumerate(zip(atom_names, rmsf)):
            rows.append(
                {
                    "mutation": r["mutation"],
                    "replicate": int(r["replicate"]),
                    "atom_name": name,
                    "atom_index": i,
                    "rmsf_angstrom": float(val),
                }
            )
    return pd.DataFrame(rows)


def plot_rmsf(rmsf: pd.DataFrame, out: Path) -> None:
    g = (
        rmsf.groupby(["mutation", "atom_name"], as_index=False)["rmsf_angstrom"]
        .agg(rmsf_mean="mean", rmsf_sem=lambda s: float(s.std(ddof=1) / np.sqrt(len(s))) if len(s) > 1 else 0.0)
    )
    g.to_csv(OUT / "tables" / "ligand_rmsf_genotype_mean_fixed.csv", index=False)
    focus = [m for m in ("WT", "Y188L", "V106I+F227C", "V106A", "G190E") if m in set(g["mutation"])]
    atoms = list(dict.fromkeys(g["atom_name"].tolist()))
    fig, ax = plt.subplots(figsize=(10, 4.2))
    x = np.arange(len(atoms))
    for mut in focus:
        sub = g[g["mutation"] == mut].set_index("atom_name").reindex(atoms)
        ax.plot(x, sub["rmsf_mean"], marker="o", ms=3, lw=1.5, label=mut)
    ax.set_xticks(x)
    ax.set_xticklabels(atoms, rotation=90, fontsize=7)
    ax.set_ylabel("DOR RMSF (Å)")
    ax.set_ylim(0, max(2.5, float(g["rmsf_mean"].max()) * 1.15))
    ax.set_title("Ligand RMSF after NNIBP Cα alignment (fixed; manual from coords)")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(alpha=0.25, linestyle=":")
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)

    # Alternative: bar of mean heavy-atom RMSF per genotype
    per = rmsf.groupby(["mutation", "replicate"], as_index=False)["rmsf_angstrom"].mean()
    summary = per.groupby("mutation", as_index=False)["rmsf_angstrom"].agg(
        rmsf_mean="mean", rmsf_sem=lambda s: float(s.std(ddof=1) / np.sqrt(len(s))) if len(s) > 1 else 0.0
    )
    summary = summary.sort_values("mutation", key=lambda s: s.map(_mutation_sort_key))
    fig, ax = plt.subplots(figsize=(10, 4))
    x = np.arange(len(summary))
    ax.bar(x, summary["rmsf_mean"], yerr=summary["rmsf_sem"], capsize=2, color="#2c6fbb", ecolor="0.4")
    ax.set_xticks(x)
    ax.set_xticklabels(summary["mutation"], rotation=90, fontsize=8)
    ax.set_ylabel("Mean DOR heavy-atom RMSF (Å)")
    ax.set_title("Ligand flexibility summary (mean ± SEM over atoms, then reps)")
    ax.grid(alpha=0.25, axis="y", linestyle=":")
    fig.tight_layout()
    fig.savefig(out.with_name("ligand_rmsf_mean_by_genotype.png"), dpi=200)
    plt.close(fig)


def plot_hbonds(hb_delta: pd.DataFrame, hb_rep: pd.DataFrame, out_dir: Path) -> None:
    # Mask self-mutation sites
    masked = hb_delta.copy()
    drop_rows = []
    for i, row in masked.iterrows():
        sites = MUTATED_SITE.get(row["mutation"], set())
        if int(row["auth_resid"]) in sites:
            drop_rows.append(i)
    masked_plot = masked.drop(index=drop_rows)

    keep_res = sorted(set(masked_plot.loc[masked_plot["occupancy_mean"] > 0.05, "auth_resid"].astype(int)))
    muts = sorted([m for m in masked_plot["mutation"].unique() if m != "WT"], key=_mutation_sort_key)
    mat = np.full((len(muts), len(keep_res)), np.nan)
    for i, mut in enumerate(muts):
        sub = masked_plot[masked_plot["mutation"] == mut].set_index("auth_resid")
        for j, r in enumerate(keep_res):
            if r in sub.index:
                mat[i, j] = float(sub.loc[r, "delta_vs_wt"])
    fig, ax = plt.subplots(figsize=(max(8, 0.45 * len(keep_res)), max(5, 0.35 * len(muts))))
    lim = float(np.nanmax(np.abs(mat))) if np.isfinite(mat).any() else 0.5
    lim = max(lim, 0.1)
    im = ax.imshow(mat, aspect="auto", cmap="coolwarm", vmin=-lim, vmax=lim)
    ax.set_xticks(range(len(keep_res)))
    ax.set_xticklabels([str(r) for r in keep_res], fontsize=8)
    ax.set_yticks(range(len(muts)))
    ax.set_yticklabels(muts, fontsize=8)
    ax.set_xlabel("Protein residue (auth)")
    ax.set_title("Δ DOR–protein H-bond occupancy vs WT\n(mutation-site residues excluded — e.g. no 227 for F227C)")
    fig.colorbar(im, ax=ax, fraction=0.03, label="Δ occupancy")
    fig.tight_layout()
    fig.savefig(out_dir / "dor_hbond_delta_heatmap.png", dpi=200)
    plt.close(fig)

    # K103 confirmation
    k = hb_delta[hb_delta["auth_resid"] == 103].copy()
    k = k.sort_values("occupancy_mean")
    fig, ax = plt.subplots(figsize=(9, 4.2))
    colors = ["#c0392b" if m in ("F227C", "A98G+F227C", "V106I+F227C") else "#2c6fbb" for m in k["mutation"]]
    ax.barh(k["mutation"], k["occupancy_mean"], xerr=k["occupancy_sem"], color=colors, ecolor="0.4", capsize=2)
    ax.axvline(float(k.loc[k["mutation"] == "WT", "occupancy_mean"].iloc[0]), color="0.2", ls="--", lw=1.2, label="WT")
    ax.set_xlabel("Lys103–DOR H-bond occupancy")
    ax.set_title("Load-bearing Lys103 H-bond (F227C family in red — not the big losers)")
    ax.set_xlim(0, 1.05)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "dor_hbond_K103_occupancy.png", dpi=200)
    plt.close(fig)

    # Residue 227 occupancy — shows the F227C drop is the mutated site itself
    r227 = hb_delta[hb_delta["auth_resid"] == 227].sort_values("occupancy_mean")
    if not r227.empty:
        fig, ax = plt.subplots(figsize=(9, 4.2))
        colors = ["#c0392b" if "F227C" in m else "#2c6fbb" for m in r227["mutation"]]
        ax.barh(r227["mutation"], r227["occupancy_mean"], color=colors)
        ax.set_xlabel("Residue 227–DOR H-bond occupancy")
        ax.set_title("Why F227C looks blue at column 227: the partner residue was mutated (Cys ≠ Phe)")
        fig.tight_layout()
        fig.savefig(out_dir / "dor_hbond_res227_occupancy.png", dpi=200)
        plt.close(fig)

    # Total pocket H-bond occupancy sum excluding mutated sites
    rows = []
    for mut, sub in hb_rep.groupby("mutation"):
        sites = MUTATED_SITE.get(mut, set())
        use = sub[~sub["auth_resid"].isin(sites)]
        # sum occupancy per rep then mean
        per = use.groupby("replicate")["occupancy"].sum()
        rows.append({"mutation": mut, "sum_occ_mean": float(per.mean()), "sum_occ_sem": float(per.std(ddof=1) / np.sqrt(len(per))) if len(per) > 1 else 0.0})
    tot = pd.DataFrame(rows).sort_values("sum_occ_mean")
    wt_val = float(tot.loc[tot["mutation"] == "WT", "sum_occ_mean"].iloc[0])
    tot["delta"] = tot["sum_occ_mean"] - wt_val
    tot.to_csv(OUT / "tables" / "dor_hbond_sum_occupancy_ex_mutsite.csv", index=False)
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.barh(tot["mutation"], tot["delta"], color=["#c0392b" if v < 0 else "#2c6fbb" for v in tot["delta"]])
    ax.axvline(0, color="0.3", lw=0.8)
    ax.set_xlabel("Σ residue H-bond occupancy − WT (mutation sites excluded)")
    ax.set_title("Overall DOR–protein H-bond budget vs WT")
    fig.tight_layout()
    fig.savefig(out_dir / "dor_hbond_sum_delta_vs_wt.png", dpi=200)
    plt.close(fig)


def plot_pocket_vs_fold(pocket_g: pd.DataFrame, fold: dict[str, float], out: Path) -> None:
    rows = []
    for _, r in pocket_g.iterrows():
        mut = r["mutation"]
        if mut == "WT" or mut not in fold:
            continue
        rows.append({**r.to_dict(), "fold": fold[mut], "log10_fold": math.log10(fold[mut])})
    df = pd.DataFrame(rows)
    if df.empty:
        return
    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    ax.errorbar(df["log10_fold"], df["volume_mean"], yerr=df["volume_sem"], fmt="o", capsize=3,
                color="#2c6fbb", ecolor="#9bbce0", zorder=3)
    for _, r in df.iterrows():
        ax.annotate(r["mutation"], (r["log10_fold"], r["volume_mean"]), textcoords="offset points",
                    xytext=(4, 4), fontsize=7)
    if len(df) >= 3:
        from scipy.stats import pearsonr, spearmanr
        x = df["log10_fold"].to_numpy()
        y = df["volume_mean"].to_numpy()
        r_p, p_p = pearsonr(x, y)
        rho, p_s = spearmanr(x, y)
        title = (f"NNIBP pocket volume vs experiment (n = {len(df)})\n"
                 f"Pearson R² = {r_p**2:.2f}, p = {p_p:.2g}  ·  Spearman ρ = {rho:.2f}  (no fit line)")
    else:
        title = "NNIBP pocket volume vs experiment"
    ax.set_title(title, fontsize=10)
    ax.set_xlabel(r"$\log_{10}$(experimental DOR fold reduction)")
    ax.set_ylabel("Pocket volume proxy (Å³)")
    ax.grid(alpha=0.25, linestyle=":")
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)
    df.to_csv(OUT / "tables" / "pocket_volume_vs_fold.csv", index=False)


def _kabsch(P: np.ndarray, Q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rigid transform mapping P → Q. Returns (R, t) with Q ≈ P @ R.T + t."""
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    H = Pc.T @ Qc
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    t = Q.mean(axis=0) - (P.mean(axis=0) @ R.T)
    return R, t


def build_dor_com_table(npy: Path, inv: pd.DataFrame) -> pd.DataFrame:
    """Ligand COM in a common WT NNIBP Cα frame (Kabsch from each frame's NNIBP)."""
    wt = inv[(inv["mutation"] == "WT") & (inv["replicate"] == 1)].iloc[0]
    ref = np.load(npy / f"{wt['npy_stem']}_nnibp_ca_xyz.npy")[0]  # (15, 3)
    rows = []
    for _, r in inv.iterrows():
        lig = np.load(npy / f"{r['npy_stem']}_lig_heavy_xyz.npy")
        nn = np.load(npy / f"{r['npy_stem']}_nnibp_ca_xyz.npy")
        for fi in range(len(lig)):
            R, t = _kabsch(nn[fi], ref)
            lig_com = lig[fi].mean(axis=0)
            nn_com = nn[fi].mean(axis=0)
            lig_g = lig_com @ R.T + t
            nn_g = nn_com @ R.T + t
            rel = lig_g - nn_g
            rows.append(
                {
                    "mutation": r["mutation"],
                    "replicate": int(r["replicate"]),
                    "frame_i": fi,
                    "com_x": float(lig_g[0]),
                    "com_y": float(lig_g[1]),
                    "com_z": float(lig_g[2]),
                    "rel_x": float(rel[0]),
                    "rel_y": float(rel[1]),
                    "rel_z": float(rel[2]),
                    "rel_r": float(np.linalg.norm(rel)),
                }
            )
    return pd.DataFrame(rows)


def enrich_dor_com_spherical(com: pd.DataFrame) -> pd.DataFrame:
    """Add spherical coords; polar axis = WT mean (lig−NNIBP) COM direction.

    alpha_deg — angle from WT mean COM axis (0 = canonical pocket location)
    beta_deg  — azimuth around that axis (−180…180)
    theta_deg / phi_deg — ordinary spherical in Kabsch WT frame (polar from +z)
    """
    out = com.copy()
    v = out[["rel_x", "rel_y", "rel_z"]].to_numpy(dtype=float)
    r = out["rel_r"].to_numpy(dtype=float)
    r_safe = np.maximum(r, 1e-8)
    out["theta_deg"] = np.degrees(np.arccos(np.clip(v[:, 2] / r_safe, -1.0, 1.0)))
    out["phi_deg"] = np.degrees(np.arctan2(v[:, 1], v[:, 0]))

    wt = out[out["mutation"] == "WT"]
    if wt.empty:
        raise ValueError("need WT frames to define spherical polar axis")
    mu = wt[["rel_x", "rel_y", "rel_z"]].mean().to_numpy(dtype=float)
    mu = mu / np.linalg.norm(mu)
    vhat = v / r_safe[:, None]
    out["alpha_deg"] = np.degrees(np.arccos(np.clip(vhat @ mu, -1.0, 1.0)))

    # Orthonormal basis with ez = μ
    ez = mu
    tmp = np.array([0.0, 0.0, 1.0]) if abs(ez[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    ex = np.cross(tmp, ez)
    ex /= np.linalg.norm(ex)
    ey = np.cross(ez, ex)
    out["beta_deg"] = np.degrees(np.arctan2(v @ ey, v @ ex))
    out.attrs["wt_com_axis"] = mu
    return out


def plot_dor_com_histograms(com: pd.DataFrame, out_dir: Path) -> None:
    """Per-genotype histograms of ligand COM (pocket-relative + absolute in WT frame)."""
    muts = sorted(com["mutation"].unique(), key=_mutation_sort_key)
    # Shared axes: relative COM components (Å) vs NNIBP Cα COM
    comps = [("rel_x", "Δx"), ("rel_y", "Δy"), ("rel_z", "Δz")]
    pad = 0.15
    lims = {
        c: (float(com[c].min()) - pad, float(com[c].max()) + pad) for c, _ in comps
    }
    r_lim = (float(com["rel_r"].min()) - pad, float(com["rel_r"].max()) + pad)
    wt_med = com.loc[com["mutation"] == "WT", "rel_r"].median() if "WT" in set(com["mutation"]) else None
    colors = {"rel_x": "#1f77b4", "rel_y": "#2ca02c", "rel_z": "#d62728"}

    def _grid(mut_list: list[str], out_path: Path, ncols: int) -> None:
        nrows = int(np.ceil(len(mut_list) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(3.5 * ncols, 3.0 * nrows), squeeze=False)
        for i, mut in enumerate(mut_list):
            ax = axes[i // ncols][i % ncols]
            sub = com[com["mutation"] == mut]
            for col, lab in comps:
                ax.hist(
                    sub[col],
                    bins=28,
                    range=lims[col],
                    density=True,
                    histtype="step",
                    linewidth=1.6,
                    color=colors[col],
                    label=lab,
                )
            ax.set_title(mut, fontweight="bold", fontsize=10)
            ax.set_xlabel("ligand COM − NNIBP COM (Å)")
            ax.set_ylabel("density")
            ax.grid(alpha=0.25, linestyle=":")
            if i == 0:
                ax.legend(fontsize=7, frameon=False, loc="upper right")
        for j in range(len(mut_list), nrows * ncols):
            axes[j // ncols][j % ncols].axis("off")
        fig.suptitle(
            "DOR ligand COM in common WT NNIBP frame\n"
            "(Kabsch on NNIBP Cα each frame; Δ = lig COM − pocket COM)",
            fontweight="bold",
        )
        fig.tight_layout()
        fig.savefig(out_path, dpi=200)
        plt.close(fig)

    def _radial_grid(mut_list: list[str], out_path: Path, ncols: int) -> None:
        nrows = int(np.ceil(len(mut_list) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(3.5 * ncols, 3.0 * nrows), squeeze=False)
        for i, mut in enumerate(mut_list):
            ax = axes[i // ncols][i % ncols]
            sub = com[com["mutation"] == mut]
            ax.hist(sub["rel_r"], bins=28, range=r_lim, density=True, color="#4c72b0", alpha=0.75, edgecolor="white", linewidth=0.4)
            if wt_med is not None:
                ax.axvline(wt_med, color="k", ls="--", lw=1.0, label="WT median" if i == 0 else None)
            ax.set_title(mut, fontweight="bold", fontsize=10)
            ax.set_xlabel("|lig COM − NNIBP COM| (Å)")
            ax.set_ylabel("density")
            ax.set_xlim(*r_lim)
            ax.grid(alpha=0.25, linestyle=":")
            if i == 0 and wt_med is not None:
                ax.legend(fontsize=7, frameon=False)
        for j in range(len(mut_list), nrows * ncols):
            axes[j // ncols][j % ncols].axis("off")
        fig.suptitle(
            "DOR binding depth proxy: |ligand COM − NNIBP Cα COM|\n"
            "(same WT-frame Kabsch; dashed = WT median)",
            fontweight="bold",
        )
        fig.tight_layout()
        fig.savefig(out_path, dpi=200)
        plt.close(fig)

    story = [m for m in STORY if m in set(com["mutation"])]
    _grid(muts, out_dir / "dor_com_hist_by_genotype.png", ncols=4)
    _grid(story, out_dir / "dor_com_hist_stories.png", ncols=3)
    _radial_grid(muts, out_dir / "dor_com_radial_hist_by_genotype.png", ncols=4)
    _radial_grid(story, out_dir / "dor_com_radial_hist_stories.png", ncols=3)

    # Summary: mean ± SEM of radial depth across replicates
    rep = (
        com.groupby(["mutation", "replicate"], as_index=False)
        .agg(rel_r_mean=("rel_r", "mean"), rel_x=("rel_x", "mean"), rel_y=("rel_y", "mean"), rel_z=("rel_z", "mean"))
    )
    summ = (
        rep.groupby("mutation", as_index=False)
        .agg(
            rel_r=("rel_r_mean", "mean"),
            rel_r_sem=("rel_r_mean", lambda s: float(s.std(ddof=1) / np.sqrt(len(s))) if len(s) > 1 else 0.0),
            rel_x=("rel_x", "mean"),
            rel_y=("rel_y", "mean"),
            rel_z=("rel_z", "mean"),
        )
    )
    summ = summ.sort_values("mutation", key=lambda s: s.map(_mutation_sort_key))
    summ.to_csv(OUT / "tables" / "dor_com_genotype_summary.csv", index=False)

    fig, ax = plt.subplots(figsize=(max(8, 0.42 * len(summ)), 4.2))
    x = np.arange(len(summ))
    ax.bar(x, summ["rel_r"], yerr=summ["rel_r_sem"], color="#4c72b0", alpha=0.85, ecolor="k", capsize=2)
    if wt_med is not None:
        ax.axhline(float(summ.loc[summ["mutation"] == "WT", "rel_r"].iloc[0]), color="k", ls="--", lw=1.0)
    ax.set_xticks(x)
    ax.set_xticklabels(summ["mutation"], rotation=55, ha="right", fontsize=8)
    ax.set_ylabel("|lig COM − NNIBP COM| (Å)")
    ax.set_title("DOR pocket depth by genotype (mean ± SEM over replicates)")
    ax.grid(axis="y", alpha=0.25, linestyle=":")
    fig.tight_layout()
    fig.savefig(out_dir / "dor_com_radial_by_genotype.png", dpi=200)
    plt.close(fig)

    plot_dor_com_spherical_2d(com, out_dir)


def plot_dor_com_spherical_2d(com: pd.DataFrame, out_dir: Path) -> None:
    """2D histograms: depth r vs polar angle α (and azimuth β) in WT-COM spherical frame."""
    if "alpha_deg" not in com.columns:
        com = enrich_dor_com_spherical(com)
    muts = sorted(com["mutation"].unique(), key=_mutation_sort_key)
    story = [m for m in STORY if m in set(com["mutation"])]

    r_pad = 0.1
    r_lim = (float(com["rel_r"].min()) - r_pad, float(com["rel_r"].max()) + r_pad)
    # Cap alpha extent at 95th percentile so outliers don't squash the map
    a_hi = float(np.percentile(com["alpha_deg"], 99))
    a_lim = (0.0, max(40.0, a_hi))
    b_lim = (-180.0, 180.0)

    wt = com[com["mutation"] == "WT"]
    wt_r = float(wt["rel_r"].median()) if len(wt) else None
    wt_a = float(wt["alpha_deg"].median()) if len(wt) else None

    def _hex_grid(
        mut_list: list[str],
        *,
        ycol: str,
        ylim: tuple[float, float],
        ylabel: str,
        out_path: Path,
        ncols: int,
        title: str,
    ) -> None:
        nrows = int(np.ceil(len(mut_list) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(3.6 * ncols, 3.2 * nrows), squeeze=False)
        for i, mut in enumerate(mut_list):
            ax = axes[i // ncols][i % ncols]
            sub = com[com["mutation"] == mut]
            hb = ax.hexbin(
                sub["rel_r"],
                sub[ycol],
                gridsize=28,
                cmap="viridis",
                mincnt=1,
                extent=(r_lim[0], r_lim[1], ylim[0], ylim[1]),
            )
            if wt_r is not None and wt_a is not None and ycol == "alpha_deg":
                ax.axvline(wt_r, color="w", ls="--", lw=0.8, alpha=0.7)
                ax.axhline(wt_a, color="w", ls="--", lw=0.8, alpha=0.7)
            ax.set_xlim(*r_lim)
            ax.set_ylim(*ylim)
            ax.set_title(mut, fontweight="bold", fontsize=10)
            ax.set_xlabel(r"$r$ (Å)")
            ax.set_ylabel(ylabel)
            fig.colorbar(hb, ax=ax, fraction=0.046, label="frames")
        for j in range(len(mut_list), nrows * ncols):
            axes[j // ncols][j % ncols].axis("off")
        fig.suptitle(title, fontweight="bold")
        fig.tight_layout()
        fig.savefig(out_path, dpi=200)
        plt.close(fig)

    title_a = (
        r"DOR COM spherical map: depth $r$ vs polar angle $\alpha$" + "\n"
        r"($\alpha$ = angle from WT mean COM axis; white dashed = WT medians)"
    )
    title_b = (
        r"DOR COM spherical map: depth $r$ vs azimuth $\beta$" + "\n"
        r"($\beta$ around WT mean COM axis; same $r$ axis across panels)"
    )
    _hex_grid(
        muts,
        ycol="alpha_deg",
        ylim=a_lim,
        ylabel=r"$\alpha$ (deg from WT axis)",
        out_path=out_dir / "dor_com_r_alpha_by_genotype.png",
        ncols=4,
        title=title_a,
    )
    _hex_grid(
        story,
        ycol="alpha_deg",
        ylim=a_lim,
        ylabel=r"$\alpha$ (deg from WT axis)",
        out_path=out_dir / "dor_com_r_alpha_stories.png",
        ncols=3,
        title=title_a,
    )
    _hex_grid(
        muts,
        ycol="beta_deg",
        ylim=b_lim,
        ylabel=r"$\beta$ (deg azimuth)",
        out_path=out_dir / "dor_com_r_beta_by_genotype.png",
        ncols=4,
        title=title_b,
    )
    _hex_grid(
        story,
        ycol="beta_deg",
        ylim=b_lim,
        ylabel=r"$\beta$ (deg azimuth)",
        out_path=out_dir / "dor_com_r_beta_stories.png",
        ncols=3,
        title=title_b,
    )


def plot_dor_conformational_maps(npy: Path, inv: pd.DataFrame, out_dir: Path) -> None:
    """Frame-level PCA of pocket-aligned DOR coords → 2D density per genotype."""
    feats = []
    meta = []
    for _, r in inv.iterrows():
        lig = np.load(npy / f"{r['npy_stem']}_lig_heavy_xyz.npy")  # (F, L, 3)
        # center each frame on ligand COM (pose in pocket frame already)
        com = lig.mean(axis=1, keepdims=True)
        centered = (lig - com).reshape(lig.shape[0], -1)
        # subsample
        pick = np.linspace(0, len(centered) - 1, num=min(80, len(centered)), dtype=int)
        for fi in pick:
            feats.append(centered[fi])
            meta.append({"mutation": r["mutation"], "replicate": int(r["replicate"])})
    X = np.vstack(feats)
    X = X - X.mean(axis=0, keepdims=True)
    u, s, vt = np.linalg.svd(X, full_matrices=False)
    scores = u[:, :2] * s[:2]
    var = (s**2) / max(1, X.shape[0] - 1)
    ratio = var / var.sum()
    meta_df = pd.DataFrame(meta)
    meta_df["PC1"] = scores[:, 0]
    meta_df["PC2"] = scores[:, 1]
    meta_df.to_csv(OUT / "tables" / "dor_pose_pca_frames.csv", index=False)

    # Global scatter — story subset (readable legend) + full panel file
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    for mut in sorted(meta_df["mutation"].unique(), key=_mutation_sort_key):
        if mut not in STORY:
            continue
        sub = meta_df[meta_df["mutation"] == mut]
        ax.scatter(sub["PC1"], sub["PC2"], s=8, alpha=0.45, label=mut)
    ax.set_xlabel(f"DOR pose PC1 ({100*ratio[0]:.1f}%)")
    ax.set_ylabel(f"DOR pose PC2 ({100*ratio[1]:.1f}%)")
    ax.set_title("DOR conformational map (pocket-aligned, frame-level PCA)")
    ax.legend(fontsize=7, markerscale=2, ncol=2)
    ax.grid(alpha=0.25, linestyle=":")
    fig.tight_layout()
    fig.savefig(out_dir / "dor_pose_pca_scatter_stories.png", dpi=200)
    plt.close(fig)

    pad = 0.5
    xlim = (meta_df["PC1"].min() - pad, meta_df["PC1"].max() + pad)
    ylim = (meta_df["PC2"].min() - pad, meta_df["PC2"].max() + pad)

    def _density_grid(muts: list[str], out_path: Path, *, ncols: int = 4) -> None:
        nrows = int(np.ceil(len(muts) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(3.6 * ncols, 3.3 * nrows), squeeze=False)
        for i, mut in enumerate(muts):
            ax = axes[i // ncols][i % ncols]
            sub = meta_df[meta_df["mutation"] == mut]
            hb = ax.hexbin(
                sub["PC1"], sub["PC2"], gridsize=30, cmap="viridis", mincnt=1,
                extent=(xlim[0], xlim[1], ylim[0], ylim[1]),
            )
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)
            ax.set_title(mut, fontweight="bold", fontsize=10)
            ax.set_xlabel("PC1")
            ax.set_ylabel("PC2")
            fig.colorbar(hb, ax=ax, fraction=0.046, label="frames")
        for j in range(len(muts), nrows * ncols):
            axes[j // ncols][j % ncols].axis("off")
        fig.suptitle(
            f"DOR pose density in PC1–PC2 (same axes for all panels)  ·  "
            f"var {100*ratio[0]:.0f}% / {100*ratio[1]:.0f}%",
            fontweight="bold",
        )
        fig.tight_layout()
        fig.savefig(out_path, dpi=200)
        plt.close(fig)

    all_muts = sorted(meta_df["mutation"].unique(), key=_mutation_sort_key)
    story_muts = [m for m in STORY if m in set(meta_df["mutation"])]
    _density_grid(all_muts, out_dir / "dor_pose_pca_density_by_genotype.png", ncols=4)
    _density_grid(story_muts, out_dir / "dor_pose_pca_density_stories.png", ncols=3)


def plot_contact_matrices(net: pd.DataFrame, out_dir: Path) -> None:
    """Full 15×15 Δcontact heatmaps for story mutants — readable vs top-8 bars."""
    focus = [m for m in ("V106A", "V106I+F227C", "V106A+F227L", "G190E", "Y188L", "F227C") if m in set(net["mutation"])]
    labels = list(NNIBP_AUTH)
    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    axes = axes.ravel()
    for ax, mut in zip(axes, focus):
        mat = np.zeros((len(labels), len(labels)))
        sub = net[net["mutation"] == mut]
        idx = {a: i for i, a in enumerate(labels)}
        for _, r in sub.iterrows():
            i, j = idx[int(r["auth_i"])], idx[int(r["auth_j"])]
            mat[i, j] = mat[j, i] = float(r["delta_freq"])
        lim = max(0.05, float(np.max(np.abs(mat))))
        im = ax.imshow(mat, cmap="coolwarm", vmin=-lim, vmax=lim)
        ax.set_xticks(range(len(labels)))
        ax.set_yticks(range(len(labels)))
        ax.set_xticklabels(labels, fontsize=6, rotation=90)
        ax.set_yticklabels(labels, fontsize=6)
        ax.set_title(f"{mut} − WT", fontsize=10, fontweight="bold")
        fig.colorbar(im, ax=ax, fraction=0.046)
    for ax in axes[len(focus):]:
        ax.axis("off")
    fig.suptitle("NNIBP residue–residue contact-frequency change vs WT\n(red = more often in contact; blue = less)", fontweight="bold")
    fig.tight_layout()
    fig.savefig(out_dir / "nnibp_contact_delta_matrices.png", dpi=200)
    plt.close(fig)


def plot_dccm_explainer(npy: Path, out: Path) -> None:
    """Side-by-side: WT DCCM, V106A DCCM, Δ — with a text caption on the figure."""
    wt = npy / "wt_mean_nnibp_dccm.npy"
    dV = npy / "dccm_delta_V106A_minus_wt.npy"
    if not (wt.is_file() and dV.is_file()):
        return
    # reconstruct V106A mean from delta
    wt_m = np.load(wt)
    delta = np.load(dV)
    mut_m = wt_m + delta
    labels = list(NNIBP_AUTH)
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.2))
    for ax, mat, title, cmap, vmin, vmax in (
        (axes[0], wt_m, "WT DCCM", "RdBu_r", -1, 1),
        (axes[1], mut_m, "V106A DCCM", "RdBu_r", -1, 1),
        (axes[2], delta, "Δ (V106A − WT)", "coolwarm", -float(np.max(np.abs(delta))), float(np.max(np.abs(delta)))),
    ):
        im = ax.imshow(mat, cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_xticks(range(len(labels)))
        ax.set_yticks(range(len(labels)))
        ax.set_xticklabels(labels, fontsize=6, rotation=90)
        ax.set_yticklabels(labels, fontsize=6)
        ax.set_title(title, fontweight="bold")
        fig.colorbar(im, ax=ax, fraction=0.046)
    fig.suptitle(
        "DCCM = pairwise correlation of Cα motions  (+1 together, −1 opposite)\n"
        "Δ heatmap = how those correlations change in the mutant — not a contact map",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def plot_nnibp_motion_pca(npy: Path, inv: pd.DataFrame, out: Path) -> None:
    """Frame-level NNIBP Cα PCA — axes of pocket motion, density by genotype."""
    feats, meta = [], []
    for _, r in inv.iterrows():
        xyz = np.load(npy / f"{r['npy_stem']}_nnibp_ca_xyz.npy")  # (F, 15, 3)
        mean = xyz.mean(axis=0)
        # align each frame roughly by subtracting mean structure (already pocket-ish)
        centered = (xyz - mean).reshape(xyz.shape[0], -1)
        pick = np.linspace(0, len(centered) - 1, num=min(60, len(centered)), dtype=int)
        for fi in pick:
            feats.append(centered[fi])
            meta.append(r["mutation"])
    X = np.vstack(feats)
    X -= X.mean(0)
    u, s, vt = np.linalg.svd(X, full_matrices=False)
    scores = u[:, :2] * s[:2]
    ratio = (s**2) / (s**2).sum()
    df = pd.DataFrame({"mutation": meta, "PC1": scores[:, 0], "PC2": scores[:, 1]})
    df.to_csv(OUT / "tables" / "nnibp_motion_pca_frames.csv", index=False)

    muts = [m for m in STORY if m in set(df["mutation"])]
    ncols = 3
    nrows = int(np.ceil(len(muts) / ncols))
    xlim = (df.PC1.min() - 0.5, df.PC1.max() + 0.5)
    ylim = (df.PC2.min() - 0.5, df.PC2.max() + 0.5)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.0 * ncols, 3.6 * nrows), squeeze=False)
    for i, mut in enumerate(muts):
        ax = axes[i // ncols][i % ncols]
        sub = df[df["mutation"] == mut]
        hb = ax.hexbin(sub.PC1, sub.PC2, gridsize=28, cmap="viridis", mincnt=1,
                       extent=(xlim[0], xlim[1], ylim[0], ylim[1]))
        ax.set_xlim(*xlim); ax.set_ylim(*ylim)
        ax.set_title(mut, fontweight="bold")
        ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
        fig.colorbar(hb, ax=ax, fraction=0.046)
    for j in range(len(muts), nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")
    fig.suptitle(
        f"NNIBP Cα motion map (frame PCA)  ·  PC1 {100*ratio[0]:.0f}% / PC2 {100*ratio[1]:.0f}%\n"
        "Same axes across panels — where each genotype samples pocket conformational space",
        fontweight="bold",
    )
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite-dir", type=Path, default=OUT)
    args = ap.parse_args()
    out = args.suite_dir
    tables, plots, npy = out / "tables", out / "plots", out / "npy"
    inv = pd.read_csv(tables / "replicate_inventory.csv")

    print("Recomputing RMSF from npy…")
    rmsf = recompute_rmsf_from_npy(npy, inv)
    rmsf.to_csv(tables / "ligand_rmsf_per_atom_fixed.csv", index=False)
    plot_rmsf(rmsf, plots / "ligand_rmsf_by_genotype.png")
    print("  C1x mean RMSF now:", float(rmsf.loc[rmsf.atom_name == "C1x", "rmsf_angstrom"].mean()))

    print("H-bond plots…")
    hb_delta = pd.read_csv(tables / "dor_hbond_occupancy_delta_vs_wt.csv")
    hb_rep = pd.read_csv(tables / "dor_hbond_occupancy_per_rep.csv")
    plot_hbonds(hb_delta, hb_rep, plots)

    print("Pocket vs fold…")
    pocket_g = pd.read_csv(tables / "pocket_volume_genotype.csv")
    fold = load_experimental(EXPERIMENTAL_CSV)
    plot_pocket_vs_fold(pocket_g, fold, plots / "pocket_volume_vs_experiment.png")

    print("DOR pose maps…")
    plot_dor_conformational_maps(npy, inv, plots)

    print("DOR COM histograms…")
    com = enrich_dor_com_spherical(build_dor_com_table(npy, inv))
    com.to_csv(tables / "dor_com_frames.csv", index=False)
    plot_dor_com_histograms(com, plots)

    print("Contact matrices…")
    net = pd.read_csv(tables / "nnibp_contact_network_delta_vs_wt.csv")
    plot_contact_matrices(net, plots)

    print("DCCM explainer + NNIBP motion PCA…")
    plot_dccm_explainer(npy, plots / "dccm_explainer_V106A.png")
    plot_nnibp_motion_pca(npy, inv, plots / "nnibp_motion_pca_density_by_genotype.png")

    # Drop misleading old contact bar figure name by rewriting a pointer note
    (plots / "README_plots.md").write_text(
        """# How to read these plots

- `dor_hbond_delta_heatmap.png` — mutation-site residues excluded. F227C’s blue square at 227 was an artifact of mutating the partner residue; see `dor_hbond_res227_occupancy.png` + `dor_hbond_K103_occupancy.png`.
- `ligand_rmsf_by_genotype.png` — **fixed**. Old C1x spike was an `md.rmsf` bug; true CF3/core RMSF is ~0.5–1.5 Å.
- `dor_pose_pca_density_by_genotype.png` — frame-level DOR conformational map (what pose clusters were trying to say).
- `dor_com_hist_by_genotype.png` / `dor_com_radial_*` — ligand COM after Kabsch to WT NNIBP; Δ = lig COM − pocket COM (the translation PCA would have eaten).
- `nnibp_motion_pca_density_by_genotype.png` — frame-level pocket Cα motion map (replaces the old per-rep mean-structure PCA scatter).
- `nnibp_contact_delta_matrices.png` — replaces the hard-to-read top-8 bar panel.
- `dccm_explainer_V106A.png` — what DCCM / ΔDCCM means.
- `pocket_volume_vs_experiment.png` — volume vs log10(fold), no fit line.
"""
    )
    print(f"Wrote updated plots under {plots}")
    return 0


if __name__ == "__main__":
    # allow running as module or script
    import sys
    if str(REPO) not in sys.path:
        sys.path.insert(0, str(REPO))
    raise SystemExit(main())
