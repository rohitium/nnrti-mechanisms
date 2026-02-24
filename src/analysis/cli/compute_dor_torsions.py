#!/usr/bin/env python3
"""Compute doravirine (2KW) torsion angles across MD trajectories.

Four torsion angles capture DOR's conformational flexibility in the NNBP,
analogous to the τ1–τ5 analysis of rilpivirine in Das et al. (PNAS 2008).

Atom names use the static topology atom name mapping derived from the 4NCG
crystal structure (element-aware nearest-neighbour matching):

    τ1  C12x – C2x  – O1x  – C9x   (pyridinone–O ether bond)
    τ2  C2x  – O1x  – C9x  – C10x  (O–chlorocyanobenzene bond)
    τ3  C4x  – N2x  – C15x – C13x  (pyridinone-N – CH₂ linker)
    τ4  N2x  – C15x – C13x – N5x   (CH₂ linker – triazolone bond)

Outputs:
    results/dor_torsions.csv
    results/plots/dor_torsions_by_mutation.png
"""
from __future__ import annotations

import json
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm

warnings.filterwarnings("ignore")

REPO      = Path(__file__).resolve().parents[3]
OUT_CSV   = REPO / "results" / "dor_torsions.csv"
PLOTS_DIR = REPO / "results" / "plots"

LIGAND_RESNAME = "2KW"

# Static torsion definitions: (name, [atom1, atom2, atom3, atom4])
# Atom names are from the analysis topology PDB (GAFF-typed, post-minimisation).
TORSIONS: list[tuple[str, list[str]]] = [
    ("tau1", ["C12x", "C2x",  "O1x",  "C9x" ]),  # pyridinone-O ether bond
    ("tau2", ["C2x",  "O1x",  "C9x",  "C10x"]),  # O-chlorocyanobenzene bond
    ("tau3", ["C4x",  "N2x",  "C15x", "C13x"]),  # pyridinone-N to CH2 linker
    ("tau4", ["N2x",  "C15x", "C13x", "N5x" ]),  # CH2 linker to triazolone
]


# ── helpers ───────────────────────────────────────────────────────────────────
def _calc_dihedral(p0: np.ndarray, p1: np.ndarray,
                   p2: np.ndarray, p3: np.ndarray) -> float:
    """Return dihedral angle in degrees (−180, +180]."""
    b1 = p1 - p0
    b2 = p2 - p1
    b3 = p3 - p2
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    m1 = np.cross(n1, b2 / (np.linalg.norm(b2) + 1e-12))
    x  = np.dot(n1, n2)
    y  = np.dot(m1, n2)
    return float(np.degrees(np.arctan2(y, x)))


def _remap(path: Path) -> Path:
    if path.exists():
        return path
    marker = "nnrti-mechanisms/"
    s = str(path)
    if marker in s:
        mapped = REPO / s.split(marker, 1)[1]
        if mapped.exists():
            return mapped
    return path


def _infer_total_ns(json_path: Path) -> float:
    try:
        j = json.loads(json_path.read_text())
        steps = int(j.get("md_production_steps_completed") or
                    j.get("md_production_steps") or 0)
        if steps > 0:
            return steps * 2.0 / 1_000_000.0
    except Exception:
        pass
    return 100.0


def _load_fold_map() -> dict:
    mf = pd.read_csv(REPO / "results" / "md_manifest.csv")
    fold = {}
    for _, row in mf.drop_duplicates("mutation").iterrows():
        v = row.get("fold_reduction")
        try:
            fold[str(row["mutation"])] = float(v) if pd.notna(v) else 1.0
        except (TypeError, ValueError):
            fold[str(row["mutation"])] = 1.0
    return fold


# ── per-replicate processing ──────────────────────────────────────────────────
def process_replicate(row: pd.Series) -> list[dict]:
    import MDAnalysis as mda

    mut_key   = str(row["safe_label"])
    mutation  = str(row["mutation"])
    replicate = int(row["replicate"])

    out_json = _remap(Path(str(row.get("output_json", "")).strip()))
    topo     = _remap(Path(str(row.get("analysis_topology_pdb", "")).strip()))
    dcd_raw  = str(row.get("analysis_dcd", "")).strip()
    dcd      = _remap(Path(dcd_raw)) if dcd_raw else Path("")

    # .bak fallback
    if not dcd.exists():
        rep_dir = topo.parent if topo.exists() else Path(".")
        for suf in (f"{mut_key}_rep{replicate:02d}_analysis.10ns.bak",
                    f"{mut_key}_rep{replicate:02d}_analysis.dcd.bak"):
            cand = _remap(rep_dir / suf)
            if cand.exists():
                dcd = cand
                break

    if not topo.exists() or not dcd.exists():
        print(f"  MISSING: {mutation} rep{replicate}")
        return []

    total_ns = _infer_total_ns(out_json) if out_json.exists() else 100.0
    fmt_kw   = {"format": "DCD"} if str(dcd).endswith(".bak") else {}
    u = mda.Universe(str(topo), str(dcd), **fmt_kw)
    n_frames = len(u.trajectory)

    # Resolve atom selections once
    lig = u.select_atoms(f"resname {LIGAND_RESNAME}")
    try:
        atom_groups = []
        for _name, atoms in TORSIONS:
            sels = [lig.select_atoms(f"name {a}") for a in atoms]
            missing = [atoms[i] for i, s in enumerate(sels) if s.n_atoms == 0]
            if missing:
                raise ValueError(f"Atoms not found: {missing}")
            atom_groups.append(sels)
    except ValueError as exc:
        print(f"  WARNING ({mutation} rep{replicate}): {exc} — skipping")
        return []

    rows_out: list[dict] = []
    for ts in u.trajectory:
        time_ns = (float(ts.frame) / max(1, n_frames - 1)) * total_ns
        base = {
            "mutation":  mutation,
            "safe_label": mut_key,
            "replicate": replicate,
            "frame":     int(ts.frame),
            "time_ns":   time_ns,
        }
        for (tname, _), sels in zip(TORSIONS, atom_groups):
            p = [s.positions[0] for s in sels]
            base[tname] = _calc_dihedral(*p)
        rows_out.append(base)

    return rows_out


# ── main ──────────────────────────────────────────────────────────────────────
def main() -> None:
    from ..result_collector import collect_md_results
    run_df = collect_md_results(REPO / "results" / "md_manifest.csv")
    FOLD = _load_fold_map()

    all_rows: list[dict] = []
    for _, row in run_df.iterrows():
        result = process_replicate(row)
        all_rows.extend(result)
        if result:
            mut = result[0]["mutation"]
            rep = result[0]["replicate"]
            vals = {t: [r[t] for r in result] for t, _ in TORSIONS}
            print(f"  {mut} rep{rep}: " +
                  "  ".join(f"{t}={np.mean(v):.0f}±{np.std(v):.0f}°"
                             for t, v in vals.items()))

    if not all_rows:
        print("No torsion data collected.")
        return

    df = pd.DataFrame(all_rows)
    df.to_csv(OUT_CSV, index=False)
    print(f"\nWrote {OUT_CSV} ({len(df)} rows)")

    plot(df, FOLD)


# ── plotting ──────────────────────────────────────────────────────────────────
def plot(df: pd.DataFrame, FOLD: dict | None = None) -> None:
    if FOLD is None:
        FOLD = _load_fold_map()
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    mutations = sorted(
        df["mutation"].unique(),
        key=lambda m: (m != "WT", FOLD.get(m, 9999.0), m),
    )
    tau_names = [t for t, _ in TORSIONS]
    n_mut = len(mutations)

    colors = dict(zip(mutations, cm.RdYlGn_r(np.linspace(0.1, 0.9, n_mut))))

    def _fold_txt(m: str) -> str:
        v = FOLD.get(m)
        if isinstance(v, float) and np.isfinite(v):
            return f"{v:.0f}×"
        return "?"

    # ── Plot 1: smooth histogram distributions, one panel per torsion ──
    bins = np.linspace(-180, 180, 73)  # 5° bins
    fig, axes = plt.subplots(len(tau_names), 1,
                             figsize=(10, 3.0 * len(tau_names)),
                             sharex=True)
    if len(tau_names) == 1:
        axes = [axes]

    for ax_idx, (ax, tau) in enumerate(zip(axes, tau_names)):
        for mut in mutations:
            vals = df[df["mutation"] == mut][tau].dropna().to_numpy()
            if len(vals) == 0:
                continue
            counts, edges = np.histogram(vals, bins=bins, density=True)
            centers = 0.5 * (edges[:-1] + edges[1:])
            ax.plot(centers, counts, color=colors[mut], linewidth=1.4,
                    label=f"{mut} ({_fold_txt(mut)})")
            ax.fill_between(centers, counts, alpha=0.07, color=colors[mut])
        ax.set_ylabel(f"{tau} density", fontsize=8)
        ax.axvline(0, color="gray", lw=0.5, ls=":")
        ax.grid(alpha=0.2, linestyle=":")
        ax.set_xlim(-185, 185)
        if ax_idx == 0:
            ax.legend(fontsize=7, frameon=False, ncol=3, loc="upper left")

    axes[-1].set_xlabel("Torsion angle (°)", fontsize=9)

    tau_desc = "τ1=pyridinone–O  |  τ2=O–Ar  |  τ3=pyr-N–CH₂  |  τ4=CH₂–triazolone"
    fig.suptitle(f"Doravirine torsion angle distributions across resistance panel\n{tau_desc}",
                 fontsize=10, fontweight="bold")
    fig.tight_layout()
    out = PLOTS_DIR / "dor_torsions_by_mutation.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}")

    # ── Plot 2: timeseries for each mutation (τ1 + τ2 only) ──
    ncols = 5
    nrows = int(np.ceil(n_mut / ncols))
    for tau_pair, suffix in [(["tau1", "tau2"], "ether"), (["tau3", "tau4"], "linker")]:
        fig2, axes2 = plt.subplots(nrows * len(tau_pair), ncols,
                                   figsize=(3.5 * ncols, 2.5 * nrows * len(tau_pair)),
                                   sharey="row")
        axes2 = np.array(axes2).reshape(nrows * len(tau_pair), ncols)
        cmap_ts = cm.get_cmap("tab10")

        for mi, mut in enumerate(mutations):
            col = mi % ncols
            base_row = (mi // ncols) * len(tau_pair)
            sub = df[df["mutation"] == mut]
            fold = FOLD.get(mut, "?")
            title = (f"{mut} ({fold:.0f}×)" if isinstance(fold, float) else mut)

            for ti, tau in enumerate(tau_pair):
                ax = axes2[base_row + ti, col]
                for rep_idx, (_, grp) in enumerate(sub.groupby("replicate")):
                    ax.plot(grp["time_ns"], grp[tau],
                            lw=0.5, alpha=0.7, color=cmap_ts(rep_idx % 10))
                ax.set_ylim(-185, 185)
                ax.axhline(0, color="gray", lw=0.3, ls=":")
                ax.tick_params(labelsize=6)
                ax.grid(axis="y", alpha=0.2)
                if col == 0:
                    ax.set_ylabel(f"{tau} (°)", fontsize=7)
                if ti == len(tau_pair) - 1:
                    ax.set_xlabel("Time (ns)", fontsize=6)
                if ti == 0:
                    ax.set_title(title, fontsize=7, fontweight="bold")

        # hide unused panels
        for mi in range(len(mutations), nrows * ncols):
            col = mi % ncols
            base_row = (mi // ncols) * len(tau_pair)
            for ti in range(len(tau_pair)):
                axes2[base_row + ti, col].set_visible(False)

        fig2.suptitle(f"Doravirine {suffix} torsion timeseries ({', '.join(tau_pair)})",
                      fontsize=10, fontweight="bold")
        fig2.tight_layout()
        out2 = PLOTS_DIR / f"dor_torsions_{suffix}_timeseries.png"
        fig2.savefig(out2, dpi=150, bbox_inches="tight")
        plt.close(fig2)
        print(f"Wrote {out2}")

    # Summary: mean ± std per mutation per torsion
    print("\n─── Torsion angle summary (mean ± std, °) ───")
    header = f"{'Mutation':<22} {'Fold':>5}  " + "  ".join(f"{t:>14}" for t in tau_names)
    print(header)
    for mut in mutations:
        sub = df[df["mutation"] == mut]
        fold = FOLD.get(mut, "?")
        parts = []
        for tau in tau_names:
            v = sub[tau].dropna().to_numpy()
            parts.append(f"{np.mean(v):+6.0f}±{np.std(v):4.0f}°" if len(v) else "         n/a")
        fold_str = f"{fold:.0f}" if isinstance(fold, float) else str(fold)
        print(f"  {mut:<20} {fold_str:>5}  " + "  ".join(parts))


if __name__ == "__main__":
    main()
