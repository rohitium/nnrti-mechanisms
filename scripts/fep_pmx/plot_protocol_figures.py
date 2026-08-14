#!/usr/bin/env python3
"""Per-genotype FEP protocol figures (01–05) for the panel scatter.

Order: cycle → hybrid → NEQ work → λ profile → Crooks estimators.

    python3 scripts/fep_pmx/plot_protocol_figures.py
    python3 scripts/fep_pmx/plot_protocol_figures.py --targets V106A Y188L

Defaults to genotypes plotted in panel_ddg_vs_experiment.png (panel_ddg.csv
rows with an experimental fold). Writes:

    results/analysis/fep_pmx/protocol/<SAFE>/01_….png … 05_….png
    results/analysis/fep_pmx/protocol/<SAFE>/run_config.json

V106A is also mirrored to protocol_v106a/ for the existing worked-example path.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.fep_jorgensen.mutations import MANUSCRIPT_PLANS
from scripts.fep_pmx.analyze_neq import read_work_values_kcal
from scripts.fep_pmx.config import FEP_PMX_ROOT
from scripts.fep_pmx.qc_neq import overlap_coefficient

PANEL_CSV = FEP_PMX_ROOT / "panel_ddg.csv"
START_PDB = Path("results/md_runs/wt/rep_01/assets/wt_md_rep01_start.pdb")
PROTOCOL_ROOT = FEP_PMX_ROOT / "protocol"

# V106A stick coloring (auth 106 = pdb 103)
VAL_UNIQUE = {"CG1", "CG2", "HG11", "HG12", "HG13", "HG21", "HG22", "HG23"}
ALA_DUMMY = {"HV1", "HV2"}


def _sanitize(name: str) -> str:
    return name.replace("+", "_").replace("/", "_").replace(" ", "")


def panel_genotypes(panel_csv: Path = PANEL_CSV) -> list[str]:
    """Genotypes shown on panel_ddg_vs_experiment.png (have experimental fold)."""
    rows = list(csv.DictReader(panel_csv.open()))
    out = []
    for r in rows:
        fold = (r.get("dor_fold_reduction") or "").strip()
        if not fold:
            continue
        out.append(r["genotype"])
    return out


def _parse_pdb(path: Path) -> list[dict]:
    rows = []
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        resid = line[22:26].strip()
        try:
            ri = int(resid)
        except ValueError:
            continue
        rows.append(
            {
                "name": line[12:16].strip(),
                "resn": line[17:21].strip(),
                "chain": line[21],
                "resid": ri,
                "xyz": np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
            }
        )
    return rows


def _leg_dir(leg_id: str, phase: str, replicate: int) -> Path:
    return FEP_PMX_ROOT / "legs" / leg_id / phase / f"rep_{replicate:02d}"


def _residue_map(leg_id: str) -> dict | None:
    path = _leg_dir(leg_id, "holo", 1) / "residue_map.json"
    if path.is_file():
        return json.loads(path.read_text())
    return None


def _hybrid_pdb(leg_id: str) -> Path | None:
    path = _leg_dir(leg_id, "holo", 1) / "hybrid.pdb"
    return path if path.is_file() else None


def _prepare(leg_id: str) -> dict:
    return json.loads((_leg_dir(leg_id, "holo", 1) / "neq" / "neq_prepare.json").read_text())


def _lambda_csv(leg_id: str) -> Path | None:
    path = FEP_PMX_ROOT / "lambda_profiles" / f"{leg_id}.csv"
    return path if path.is_file() else None


def _load_units(leg_id: str) -> list[dict]:
    rows = []
    for phase in ("holo", "apo"):
        for rep in (1, 2, 3):
            meta_path = _leg_dir(leg_id, phase, rep) / "neq" / "analysis" / "analysis.json"
            if not meta_path.is_file():
                continue
            meta = json.loads(meta_path.read_text())
            wf = np.array(read_work_values_kcal(Path(meta["integ_fwd"])))
            wr = np.array(read_work_values_kcal(Path(meta["integ_rev"])))
            rows.append(
                {
                    **meta,
                    "leg_id": leg_id,
                    "wf": wf,
                    "wr": wr,
                    "overlap": overlap_coefficient(wf, wr),
                }
            )
    return rows


def _cycle_means(summary: dict) -> tuple[float, float, float, float]:
    """Mean summed ΔG_holo / ΔG_apo across legs, plus ΔΔG ± SEM."""
    legs = summary["legs"]
    reps = sorted({r["replicate"] for leg in legs for r in leg["per_rep"]})
    holo_reps, apo_reps = [], []
    for rep in reps:
        h = 0.0
        a = 0.0
        ok = True
        for leg in legs:
            row = next((x for x in leg["per_rep"] if x["replicate"] == rep), None)
            if row is None:
                ok = False
                break
            h += float(row["holo_dg"])
            a += float(row["apo_dg"])
        if ok:
            holo_reps.append(h)
            apo_reps.append(a)
    return (
        float(np.mean(holo_reps)),
        float(np.mean(apo_reps)),
        float(summary["ddg_bind"]),
        float(summary["sem"]),
    )


def _plot_cycle(genotype: str, summary: dict, out_path: Path) -> None:
    holo, apo, ddg, sem = _cycle_means(summary)
    n_legs = len(summary["legs"])
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 7)
    ax.axis("off")
    boxes = {
        "wt_h": (1.1, 4.6, "WT RT · DOR\n(holo, λ = 0)"),
        "mut_h": (6.3, 4.6, f"{genotype} RT · DOR\n(holo, λ = 1)"),
        "wt_a": (1.1, 1.15, "WT RT\n(apo, λ = 0)"),
        "mut_a": (6.3, 1.15, f"{genotype} RT\n(apo, λ = 1)"),
    }
    for x, y, text in boxes.values():
        ax.add_patch(
            FancyBboxPatch(
                (x, y), 2.6, 1.45, boxstyle="round,pad=0.08,rounding_size=0.15",
                facecolor="#eef3f8", edgecolor="#1f4e79", lw=1.4,
            )
        )
        ax.text(x + 1.3, y + 0.72, text, ha="center", va="center", fontsize=10)

    def arrow(a, b, label, color="#1f4e79"):
        ax.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="-|>", color=color, lw=1.6))
        mx, my = (a[0] + b[0]) / 2, (a[1] + b[1]) / 2
        ax.text(mx, my + 0.18, label, ha="center", va="bottom", fontsize=9, color=color, fontweight="bold")

    arrow((3.7, 5.35), (6.3, 5.35), rf"$\Delta G_{{\mathrm{{holo}}}}$ = {holo:.2f}")
    arrow((3.7, 1.9), (6.3, 1.9), rf"$\Delta G_{{\mathrm{{apo}}}}$ = {apo:.2f}")
    ax.annotate("", xy=(2.4, 4.6), xytext=(2.4, 2.6),
                arrowprops=dict(arrowstyle="-|>", color="0.35", lw=1.2))
    ax.text(1.05, 3.55, r"$\Delta G_{\mathrm{bind}}^{\mathrm{WT}}$", fontsize=9, color="0.35", rotation=90, va="center")
    ax.annotate("", xy=(7.6, 4.6), xytext=(7.6, 2.6),
                arrowprops=dict(arrowstyle="-|>", color="0.35", lw=1.2))
    ax.text(8.95, 3.55, rf"$\Delta G_{{\mathrm{{bind}}}}^{{\mathrm{{{genotype}}}}}$",
            fontsize=9, color="0.35", rotation=90, va="center")
    title = f"Thermodynamic cycle  ·  WT → {genotype}"
    if n_legs > 1:
        leg_ids = " + ".join(leg["leg_id"] for leg in summary["legs"])
        title += f"\n({n_legs} additive legs: {leg_ids})"
    ax.set_title(title, fontweight="bold")
    ax.text(
        5.0, 0.35,
        rf"$\Delta\Delta G_{{\mathrm{{bind}}}} = \Delta G_{{\mathrm{{holo}}}} - \Delta G_{{\mathrm{{apo}}}} = {ddg:+.2f} \pm {sem:.2f}$ kcal/mol",
        ha="center", va="center", fontsize=11,
        bbox=dict(boxstyle="round,pad=0.35", facecolor="#fff6e5", edgecolor="#c47b16"),
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _project(xyz: np.ndarray, origin: np.ndarray, xhat: np.ndarray, yhat: np.ndarray) -> np.ndarray:
    rel = xyz - origin
    return np.column_stack([rel @ xhat, rel @ yhat])


def _draw_sticks(ax, xy, pairs, *, color, lw, ls="-", alpha=1.0):
    for a, b in pairs:
        if a not in xy or b not in xy:
            continue
        ax.plot([xy[a][0], xy[b][0]], [xy[a][1], xy[b][1]], color=color, lw=lw, ls=ls, alpha=alpha, zorder=3, solid_capstyle="round")


def _draw_atoms(ax, xy, names, *, color, s, alpha=1.0, zorder=4):
    pts = [xy[n] for n in names if n in xy]
    if not pts:
        return
    pts = np.vstack(pts)
    ax.scatter(pts[:, 0], pts[:, 1], s=s, c=color, alpha=alpha, zorder=zorder, edgecolors="0.15", linewidths=0.4)


def _plot_hybrid_v106a(out_path: Path, hybrid_pdb: Path) -> None:
    hy = _parse_pdb(hybrid_pdb)
    st = _parse_pdb(START_PDB)
    mut = [a for a in hy if a["chain"] == "A" and a["resid"] == 103]
    lig = [a for a in st if a["resn"] == "2KW" and a["name"][0] != "H"]
    lookup = {a["name"]: a["xyz"] for a in mut}
    origin = lookup["CB"]
    v1 = lookup["CG1"] - origin
    v2 = lookup["CG2"] - origin
    xhat = v1 / np.linalg.norm(v1)
    ytmp = v2 - xhat * np.dot(v2, xhat)
    yhat = ytmp / np.linalg.norm(ytmp)
    xy = {a["name"]: _project(a["xyz"][None, :], origin, xhat, yhat)[0] for a in mut}
    p_lig = _project(np.vstack([a["xyz"] for a in lig]), origin, xhat, yhat)
    xyz_lig = np.vstack([a["xyz"] for a in lig])
    backbone = [("N", "CA"), ("CA", "C"), ("C", "O"), ("CA", "CB")]
    val_bonds = [("CB", "CG1"), ("CB", "CG2")]
    ala_bonds = [("CB", "HV1"), ("CB", "HV2")]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.6), sharex=True, sharey=True)
    panels = [
        {"ax": axes[0], "title": r"$\lambda = 0$", "subtitle": "WT: Valine is on, Alanine is off",
         "show_val": True, "ghost_val": False, "show_ala": False, "ghost_ala": True},
        {"ax": axes[1], "title": r"$\lambda = 1$", "subtitle": "V106A: Valine is off, Alanine is on",
         "show_val": False, "ghost_val": True, "show_ala": True, "ghost_ala": False},
    ]
    for panel in panels:
        ax = panel["ax"]
        for i in range(len(lig)):
            for j in range(i + 1, len(lig)):
                if np.linalg.norm(xyz_lig[i] - xyz_lig[j]) < 1.85:
                    ax.plot([p_lig[i, 0], p_lig[j, 0]], [p_lig[i, 1], p_lig[j, 1]], color="#9bb6d3", lw=1.0, zorder=1)
        ax.scatter(p_lig[:, 0], p_lig[:, 1], s=12, c="#9bb6d3", zorder=1)
        _draw_sticks(ax, xy, backbone, color="#333333", lw=2.2)
        _draw_atoms(ax, xy, ["N", "CA", "C", "O", "CB"], color="#4d4d4d", s=70)
        if panel["show_val"] or panel["ghost_val"]:
            alpha = 1.0 if panel["show_val"] else 0.22
            _draw_sticks(ax, xy, val_bonds, color="#2ca02c", lw=2.8, alpha=alpha)
            _draw_atoms(ax, xy, ["CG1", "CG2"], color="#2ca02c", s=160, alpha=alpha)
        if panel["show_ala"] or panel["ghost_ala"]:
            alpha = 1.0 if panel["show_ala"] else 0.22
            _draw_sticks(ax, xy, ala_bonds, color="#d62728", lw=2.0, ls="--" if panel["ghost_ala"] else "-", alpha=alpha)
            _draw_atoms(ax, xy, ["HV1", "HV2"], color="#d62728", s=90, alpha=alpha)
        for name, label in (("CB", "Cβ"), ("CG1", "Cγ1"), ("CG2", "Cγ2"), ("CA", "Cα")):
            if name in xy:
                ax.annotate(label, xy[name], textcoords="offset points", xytext=(4, 4), fontsize=8, color="0.2")
        ax.set_aspect("equal")
        xs = [xy[n][0] for n in ("CA", "CB", "CG1", "CG2", "HV1", "HV2") if n in xy]
        ys = [xy[n][1] for n in ("CA", "CB", "CG1", "CG2", "HV1", "HV2") if n in xy]
        pad = 2.2
        ax.set_xlim(min(xs) - pad, max(xs) + pad)
        ax.set_ylim(min(ys) - pad, max(ys) + pad)
        ax.annotate("doravirine", p_lig.mean(axis=0), textcoords="offset points", xytext=(6, 8),
                    fontsize=9, color="#3d6fa0", fontstyle="italic")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(panel["title"], fontweight="bold", fontsize=11)
        ax.text(0.5, -0.04, panel["subtitle"], transform=ax.transAxes, ha="center", va="top", fontsize=9, color="0.3")
        for spine in ax.spines.values():
            spine.set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _plot_hybrid_generic(leg_id: str, genotype: str, hybrid_pdb: Path, rmap: dict, out_path: Path) -> None:
    """Overlay hybrid side chain + DOR for any leg that has a local hybrid.pdb."""
    hy = _parse_pdb(hybrid_pdb)
    st = _parse_pdb(START_PDB) if START_PDB.is_file() else []
    pdb_resid = int(rmap["pdb_residue_id"])
    mut = [a for a in hy if a["chain"] == rmap.get("chain_id", "A") and a["resid"] == pdb_resid]
    if not mut:
        _plot_hybrid_schematic(genotype, summary_legs=[leg_id], out_path=out_path)
        return
    lig = [a for a in st if a["resn"] == "2KW" and a["name"][0] != "H"]
    lookup = {a["name"]: a["xyz"] for a in mut}
    origin = lookup.get("CB", lookup.get("CA", mut[0]["xyz"]))
    # Build a local plane from CB and the two farthest heavy sidechain atoms.
    heavy = [a for a in mut if a["name"][0] != "H" and a["name"] not in {"N", "CA", "C", "O"}]
    if len(heavy) >= 2:
        d = sorted(heavy, key=lambda a: -np.linalg.norm(a["xyz"] - origin))
        v1 = d[0]["xyz"] - origin
        v2 = d[1]["xyz"] - origin
    else:
        v1 = np.array([1.0, 0.0, 0.0])
        v2 = np.array([0.0, 1.0, 0.0])
    xhat = v1 / (np.linalg.norm(v1) + 1e-12)
    ytmp = v2 - xhat * np.dot(v2, xhat)
    if np.linalg.norm(ytmp) < 1e-8:
        ytmp = np.array([0.0, 1.0, 0.0])
    yhat = ytmp / np.linalg.norm(ytmp)
    xy = {a["name"]: _project(a["xyz"][None, :], origin, xhat, yhat)[0] for a in mut}
    old = rmap.get("old_residue", "?")
    new = rmap.get("new_residue", "?")
    auth = rmap.get("auth_residue_id", "?")
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.6), sharex=True, sharey=True)
    for ax, lam, sub in (
        (axes[0], 0, f"λ = 0 · {old} on (WT)"),
        (axes[1], 1, f"λ = 1 · {new} on ({genotype})"),
    ):
        if lig:
            p_lig = _project(np.vstack([a["xyz"] for a in lig]), origin, xhat, yhat)
            xyz_lig = np.vstack([a["xyz"] for a in lig])
            for i in range(len(lig)):
                for j in range(i + 1, len(lig)):
                    if np.linalg.norm(xyz_lig[i] - xyz_lig[j]) < 1.85:
                        ax.plot([p_lig[i, 0], p_lig[j, 0]], [p_lig[i, 1], p_lig[j, 1]], color="#9bb6d3", lw=1.0, zorder=1)
            ax.scatter(p_lig[:, 0], p_lig[:, 1], s=12, c="#9bb6d3", zorder=1)
            ax.annotate("doravirine", p_lig.mean(axis=0), textcoords="offset points", xytext=(6, 8),
                        fontsize=9, color="#3d6fa0", fontstyle="italic")
        # backbone
        _draw_sticks(ax, xy, [("N", "CA"), ("CA", "C"), ("C", "O"), ("CA", "CB")], color="#333333", lw=2.2)
        _draw_atoms(ax, xy, ["N", "CA", "C", "O", "CB"], color="#4d4d4d", s=70)
        # all other heavy atoms as the morphing side chain
        side = [a["name"] for a in heavy]
        for a in heavy:
            for b in heavy:
                if a["name"] >= b["name"]:
                    continue
                if np.linalg.norm(a["xyz"] - b["xyz"]) < 1.85:
                    _draw_sticks(ax, xy, [(a["name"], b["name"])], color="#2c6fbb" if lam == 0 else "#d1642f", lw=2.2)
        _draw_atoms(ax, xy, side, color="#2c6fbb" if lam == 0 else "#d1642f", s=90)
        ax.set_aspect("equal")
        pts = np.vstack(list(xy.values()))
        pad = 2.2
        ax.set_xlim(pts[:, 0].min() - pad, pts[:, 0].max() + pad)
        ax.set_ylim(pts[:, 1].min() - pad, pts[:, 1].max() + pad)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(rf"$\lambda = {lam}$", fontweight="bold", fontsize=11)
        ax.text(0.5, -0.04, f"{sub}  ·  auth {auth}", transform=ax.transAxes,
                ha="center", va="top", fontsize=9, color="0.3")
        for spine in ax.spines.values():
            spine.set_visible(False)
    fig.suptitle(f"Hybrid topology  ·  {leg_id}", fontweight="bold")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _plot_hybrid_schematic(genotype: str, summary_legs: list[str], out_path: Path) -> None:
    """Text schematic when hybrid.pdb is not local (most Sherlock-only legs)."""
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0))
    mut_lines = []
    for leg_id in summary_legs:
        rmap = _residue_map(leg_id)
        if rmap:
            mut_lines.append(
                f"{leg_id}:  {rmap.get('old_residue')}{rmap.get('auth_residue_id')} → {rmap.get('new_residue')}"
            )
        else:
            # fall back to MutationLeg.mutation string via MANUSCRIPT_PLANS
            mut_lines.append(leg_id)
    for ax, lam, state in (
        (axes[0], 0, "A-state (WT / start) on"),
        (axes[1], 1, "B-state (mutant / end) on"),
    ):
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
        ax.add_patch(FancyBboxPatch((0.08, 0.2), 0.84, 0.55, boxstyle="round,pad=0.04,rounding_size=0.08",
                                    facecolor="#eef3f8", edgecolor="#1f4e79", lw=1.4))
        ax.text(0.5, 0.55, rf"$\lambda = {lam}$", ha="center", va="center", fontsize=16, fontweight="bold")
        ax.text(0.5, 0.38, state, ha="center", va="center", fontsize=10, color="0.3")
        ax.set_title(rf"$\lambda = {lam}$", fontweight="bold")
    fig.suptitle(
        f"Hybrid topology  ·  {genotype}\n"
        + "\n".join(mut_lines)
        + "\n(hybrid.pdb not local — schematic only; sticks available for V106A / Y188L)",
        fontweight="bold",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _plot_hybrid(genotype: str, summary: dict, out_path: Path) -> dict:
    leg_ids = [leg["leg_id"] for leg in summary["legs"]]
    # Sticks only from the genotype-defining (last) leg — don't silently reuse an
    # earlier leg's hybrid (e.g. wt_to_V106A for V106A+F227L).
    last = leg_ids[-1]
    hybrid = _hybrid_pdb(last)
    rmap = _residue_map(last)
    if hybrid is not None and rmap is not None:
        if genotype == "V106A" and last == "wt_to_V106A":
            _plot_hybrid_v106a(out_path, hybrid)
            return {"mode": "sticks_v106a", "leg_id": last, "hybrid_pdb": str(hybrid)}
        _plot_hybrid_generic(last, genotype, hybrid, rmap, out_path)
        return {"mode": "sticks_generic", "leg_id": last, "hybrid_pdb": str(hybrid)}
    _plot_hybrid_schematic(genotype, leg_ids, out_path)
    return {"mode": "schematic", "leg_ids": leg_ids, "hybrid_pdb": None}


def _gauss_pdf(x: np.ndarray, mu: float, sigma: float) -> np.ndarray:
    return np.exp(-0.5 * ((x - mu) / sigma) ** 2) / (sigma * np.sqrt(2.0 * np.pi))


def _gauss_intersection(mu_f: float, sig_f: float, mu_r: float, sig_r: float) -> float:
    a = 1.0 / sig_f**2 - 1.0 / sig_r**2
    b = -2.0 * mu_f / sig_f**2 + 2.0 * mu_r / sig_r**2
    c = mu_f**2 / sig_f**2 - mu_r**2 / sig_r**2 + 2.0 * np.log(sig_f / sig_r)
    if abs(a) < 1e-12:
        return float(-c / b)
    disc = b * b - 4.0 * a * c
    roots = [(-b + np.sqrt(disc)) / (2.0 * a), (-b - np.sqrt(disc)) / (2.0 * a)]
    lo, hi = sorted((mu_f, mu_r))
    between = [r for r in roots if lo <= r <= hi]
    return float(between[0] if between else roots[int(np.argmin(np.abs(np.array(roots) - 0.5 * (mu_f + mu_r))))])


def _plot_work(genotype: str, summary: dict, units_by_leg: dict[str, list[dict]], out_path: Path) -> None:
    legs = [leg["leg_id"] for leg in summary["legs"]]
    n_legs = len(legs)
    fig, axes = plt.subplots(n_legs * 2, 3, figsize=(11.2, 3.1 * n_legs * 2), sharey=False, squeeze=False)
    title_bits = []
    for li, leg_id in enumerate(legs):
        prep = _prepare(leg_id)
        title_bits.append(f"{leg_id}: {prep['equil_ns']:g} ns → {prep['n_snapshots']} snaps → {prep['switch_ps']:.0f} ps")
        units = units_by_leg[leg_id]
        for row_i, phase in enumerate(("holo", "apo")):
            row = li * 2 + row_i
            for col, rep in enumerate((1, 2, 3)):
                ax = axes[row][col]
                u = next((x for x in units if x["phase"] == phase and x["replicate"] == rep), None)
                if u is None:
                    ax.set_visible(False)
                    continue
                ax.hist(u["wf"], bins=18, alpha=0.55, density=True, color="#2c6fbb", label=r"$W_f$")
                ax.hist(u["wr"], bins=18, alpha=0.55, density=True, color="#d1642f", label=r"$W_r$")
                ax.axvline(u["bar_dg"], color="0.15", ls="--", lw=1.2, label=f"BAR {u['bar_dg']:.2f}")
                ax.set_title(f"{leg_id}  {phase} rep{rep}  ov {u['overlap']:.2f}", fontsize=8, fontweight="bold")
                ax.grid(alpha=0.2, linestyle=":")
                if col == 0:
                    ax.set_ylabel("density")
                if row == n_legs * 2 - 1:
                    ax.set_xlabel("work (kcal/mol)")
                if row == 0 and col == 2:
                    ax.legend(fontsize=7)
    fig.suptitle(f"NEQ work  ·  {genotype}\n" + "  ·  ".join(title_bits), fontweight="bold", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _plot_lambda(genotype: str, summary: dict, out_path: Path) -> dict:
    legs = [leg["leg_id"] for leg in summary["legs"]]
    available = [(leg_id, p) for leg_id in legs if (p := _lambda_csv(leg_id)) is not None]
    missing = [leg_id for leg_id in legs if _lambda_csv(leg_id) is None]
    if missing:
        fig, ax = plt.subplots(figsize=(9.6, 3.8))
        ax.axis("off")
        have = ", ".join(lid for lid, _ in available) or "(none)"
        ax.text(
            0.5, 0.5,
            f"λ profile incomplete for {genotype}\n"
            f"have: {have}\n"
            f"missing: {', '.join(missing)}\n\n"
            f"(needs switches/*/dgdl.xvg on Sherlock)\n"
            f"python3 scripts/fep_pmx/plot_lambda_profile.py --legs {' '.join(legs)}",
            ha="center", va="center", fontsize=11,
        )
        fig.savefig(out_path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        return {"status": "placeholder", "legs": legs, "have": [lid for lid, _ in available], "missing": missing}

    n = len(available)
    fig, axes = plt.subplots(n, 2, figsize=(9.6, 3.4 * n), sharey=False, squeeze=False)
    for row, (leg_id, csv_path) in enumerate(available):
        data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
        lam = data[:, 0]
        end_label = leg_id.split("_to_")[-1].replace("_", "+")
        for col, phase, i in ((0, "holo", 1), (1, "apo", 4)):
            ax = axes[row][col]
            wf, g, diss = data[:, i], data[:, i + 1], data[:, i + 2]
            ax.fill_between(lam, g, wf, color="#f4c7a8", alpha=0.7, label="dissipation")
            ax.plot(lam, wf, color="#2c6fbb", lw=1.6, label=r"$\langle W_f(\lambda)\rangle$")
            ax.plot(lam, g, color="#2ca02c", lw=1.6, label=r"$G_f(\lambda)$")
            ax.set_title(f"{leg_id}  ·  {phase}", fontweight="bold", fontsize=10)
            ax.set_xlabel(rf"$\lambda$ (0 = start, 1 = {end_label})")
            ax.grid(alpha=0.25, linestyle=":")
            if col == 0:
                ax.set_ylabel("kcal/mol")
            ax.legend(fontsize=7)
    fig.suptitle(
        rf"Free-energy profile along $\lambda$  ·  {genotype}  ·  forward switches pooled",
        fontweight="bold",
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return {"status": "ok", "legs": [lid for lid, _ in available]}


def _plot_crooks_phase(ax, units: list[dict], phase: str) -> None:
    sub = [u for u in units if u["phase"] == phase]
    wf = np.concatenate([u["wf"] for u in sub])
    wr = np.concatenate([u["wr"] for u in sub])
    bar = float(np.mean([u["bar_dg"] for u in sub]))
    cgi_pmx = float(np.mean([u["cgi_dg"] for u in sub if u.get("cgi_dg") is not None]))
    jarz = float(np.mean([u["jarz_dg_mean"] for u in sub if u.get("jarz_dg_mean") is not None]))
    ov = overlap_coefficient(wf, wr)
    ax.hist(wf, bins=28, alpha=0.35, density=True, color="#2c6fbb", label=r"$P(W_f)$")
    ax.hist(wr, bins=28, alpha=0.35, density=True, color="#d1642f", label=r"$P(W_r)$")
    mu_f, sig_f = float(wf.mean()), float(wf.std(ddof=1))
    mu_r, sig_r = float(wr.mean()), float(wr.std(ddof=1))
    x = np.linspace(min(wf.min(), wr.min()) - 0.5, max(wf.max(), wr.max()) + 0.5, 400)
    yf = _gauss_pdf(x, mu_f, sig_f)
    yr = _gauss_pdf(x, mu_r, sig_r)
    ax.plot(x, yf, color="#2c6fbb", lw=2.0, label=r"Gaussian fit $W_f$")
    ax.plot(x, yr, color="#d1642f", lw=2.0, label=r"Gaussian fit $W_r$")
    cgi = _gauss_intersection(mu_f, sig_f, mu_r, sig_r)
    y_cross = float(_gauss_pdf(np.array([cgi]), mu_f, sig_f)[0])
    ax.plot(cgi, y_cross, "o", color="#6a3d9a", ms=7, zorder=6)
    ax.axvline(cgi, color="#6a3d9a", ls=":", lw=2.0, label=rf"CGI {cgi:.2f}")
    ax.axvline(bar, color="0.1", ls="--", lw=1.8, label=rf"BAR {bar:.2f}")
    ax.axvline(jarz, color="#2ca02c", ls="-.", lw=1.8, label=rf"Jarz {jarz:.2f}")
    ax.set_title(f"{phase}  ·  ov {ov:.2f}  ·  pmx CGI {cgi_pmx:.2f}", fontsize=9, fontweight="bold")
    ax.set_xlabel(r"work $W$ (kcal/mol)")
    ax.grid(alpha=0.2, linestyle=":")
    ax.legend(fontsize=6, loc="upper left")


def _plot_crooks(genotype: str, summary: dict, units_by_leg: dict[str, list[dict]], out_path: Path) -> None:
    legs = [leg["leg_id"] for leg in summary["legs"]]
    n = len(legs)
    fig, axes = plt.subplots(n, 2, figsize=(10.6, 4.2 * n), sharey=False, squeeze=False)
    for row, leg_id in enumerate(legs):
        for col, phase in enumerate(("holo", "apo")):
            ax = axes[row][col]
            _plot_crooks_phase(ax, units_by_leg[leg_id], phase)
            if n > 1:
                ax.set_title(f"{leg_id}  ·  " + ax.get_title(), fontsize=9, fontweight="bold")
        axes[row][0].set_ylabel("density")
    ddg = float(summary["ddg_bind"])
    sem = float(summary["sem"])
    fig.suptitle(
        rf"Estimators on Crooks work histograms  ·  {genotype}  ·  "
        rf"$\Delta\Delta G_{{\mathrm{{BAR}}}} = {ddg:+.2f} \pm {sem:.2f}$ kcal/mol",
        fontweight="bold",
    )
    fig.text(
        0.5, 0.01,
        r"CGI: Gaussian ∩ of $P(W_f)$, $P(W_r)$   ·   "
        r"BAR: Crooks $P(W_f)/P(W_r)=e^{\beta(W-\Delta G)}$   ·   "
        r"Jarzynski: $\Delta G=-kT\ln\langle e^{-\beta W_f}\rangle$ (forward only)",
        ha="center", va="bottom", fontsize=8.5, color="0.25",
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_genotype(genotype: str, out_root: Path = PROTOCOL_ROOT) -> Path:
    if genotype not in MANUSCRIPT_PLANS:
        raise KeyError(f"{genotype} not in MANUSCRIPT_PLANS")
    summary_path = FEP_PMX_ROOT / "targets" / genotype / "summary.json"
    if not summary_path.is_file():
        raise FileNotFoundError(summary_path)
    summary = json.loads(summary_path.read_text())
    safe = _sanitize(genotype)
    out = out_root / safe
    out.mkdir(parents=True, exist_ok=True)

    units_by_leg = {}
    for leg in summary["legs"]:
        units_by_leg[leg["leg_id"]] = _load_units(leg["leg_id"])
        if len(units_by_leg[leg["leg_id"]]) < 6:
            print(f"  warn {genotype}/{leg['leg_id']}: only {len(units_by_leg[leg['leg_id']])}/6 analysis units")

    _plot_cycle(genotype, summary, out / "01_thermodynamic_cycle.png")
    hybrid_meta = _plot_hybrid(genotype, summary, out / "02_hybrid_topology.png")
    _plot_work(genotype, summary, units_by_leg, out / "03_neq_work.png")
    lambda_meta = _plot_lambda(genotype, summary, out / "04_lambda_profile.png")
    _plot_crooks(genotype, summary, units_by_leg, out / "05_crooks_overlap.png")

    preps = {leg["leg_id"]: _prepare(leg["leg_id"]) for leg in summary["legs"]}
    payload = {
        "genotype": genotype,
        "leg_ids": [leg["leg_id"] for leg in summary["legs"]],
        "ddg_bind": summary["ddg_bind"],
        "sem": summary["sem"],
        "fold": summary.get("fold"),
        "panels": [
            "01_thermodynamic_cycle.png",
            "02_hybrid_topology.png",
            "03_neq_work.png",
            "04_lambda_profile.png",
            "05_crooks_overlap.png",
        ],
        "hybrid": hybrid_meta,
        "lambda_profile": lambda_meta,
        "per_leg_prepare": {
            lid: {"switch_ps": p["switch_ps"], "equil_ns": p["equil_ns"], "n_snapshots": p["n_snapshots"]}
            for lid, p in preps.items()
        },
        "notes": [
            "Numbered 01–05. Same series as the V106A worked example.",
            "02 uses local hybrid.pdb sticks when present; otherwise a schematic.",
            "04 is a placeholder unless lambda_profiles/<leg>.csv exists (needs Sherlock dgdl).",
            "Multi-leg genotypes sum ΔG_holo/ΔG_apo across additive legs for the cycle.",
        ],
    }
    (out / "run_config.json").write_text(json.dumps(payload, indent=2) + "\n")

    # Keep the historical worked-example path in sync.
    if genotype == "V106A":
        legacy = FEP_PMX_ROOT / "protocol_v106a"
        legacy.mkdir(parents=True, exist_ok=True)
        for name in payload["panels"] + ["run_config.json"]:
            src = out / name
            dst = legacy / name
            dst.write_bytes(src.read_bytes())

    return out


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Per-genotype FEP protocol figure series (01–05).")
    p.add_argument("--targets", nargs="+", default=None,
                   help="Genotypes (default: panel_ddg.csv rows with experimental fold)")
    p.add_argument("--include-no-fold", action="store_true",
                   help="Also include panel_ddg.csv rows that lack dor_fold_reduction")
    p.add_argument("--output-root", type=Path, default=PROTOCOL_ROOT)
    args = p.parse_args(argv)

    if args.targets:
        targets = args.targets
    else:
        targets = panel_genotypes()
        if args.include_no_fold:
            all_rows = [r["genotype"] for r in csv.DictReader(PANEL_CSV.open())]
            targets = list(dict.fromkeys(all_rows))

    print(f"targets ({len(targets)}): {', '.join(targets)}")
    ok, fail = 0, []
    for g in targets:
        try:
            out = plot_genotype(g, out_root=args.output_root)
            print(f"wrote {out}")
            ok += 1
        except Exception as exc:  # noqa: BLE001
            print(f"FAIL {g}: {exc}")
            fail.append((g, str(exc)))
    print(f"done: {ok} ok, {len(fail)} failed")
    return 1 if fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
