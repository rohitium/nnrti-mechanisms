#!/usr/bin/env python3
"""Compute doravirine (2KW) torsion angles across MD trajectories.

Four torsion angles capture DOR's conformational flexibility in the NNBP,
analogous to the τ1–τ5 analysis of rilpivirine in Das et al. (PNAS 2008).

Atom names use the static topology atom name mapping derived from the 4NCG
crystal structure (element-aware nearest-neighbour matching):

    τ1  C12x – C2x  – O1x  – C9x   (pyridinone–O ether bond)
    τ2  C2x  – O1x  – C9x  – C10x  (O–chlorocyanobenzene bond)
    τ3  C4x  – N2x  – C15x – C14x  (pyridinone-N – CH₂ linker)
    τ4  N2x  – C15x – C14x – N5x   (CH₂ linker – triazolone bond)

Outputs:
    results/dor_torsions.csv
    results/plots/dor_torsions_by_mutation.png
"""
from __future__ import annotations

import argparse
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
    ("tau3", ["C4x",  "N2x",  "C15x", "C14x"]),  # pyridinone-N to CH2 linker
    ("tau4", ["N2x",  "C15x", "C14x", "N5x" ]),  # CH2 linker to triazolone
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


def _infer_total_ns(
    json_path: Path | None,
    *,
    fallback_state_csv: Path | None = None,
) -> float | None:
    if json_path is not None and json_path.is_file():
        try:
            j = json.loads(json_path.read_text())
            steps = int(j.get("md_production_steps_completed") or j.get("md_production_steps") or 0)
            if steps > 0:
                return steps * 2.0 / 1_000_000.0
        except Exception:
            pass

        # Fall back to matching md_state.csv beside output JSON.
        m = re.match(r"^(.+)_rep(\d{2}).*\.json$", json_path.name)
        if m:
            state_csv = json_path.parent / f"{m.group(1)}_rep{m.group(2)}_md_state.csv"
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

    if fallback_state_csv is not None and fallback_state_csv.exists():
        try:
            sdf = pd.read_csv(fallback_state_csv)
            for col in ('#"Step"', "Step"):
                if col in sdf.columns:
                    steps = pd.to_numeric(sdf[col], errors="coerce").dropna()
                    if not steps.empty:
                        return float(steps.max()) * 2.0 / 1_000_000.0
        except Exception:
            pass

    return None


def _infer_completed_steps(
    json_path: Path | None,
    *,
    fallback_state_csv: Path | None = None,
    source: str = "max",
) -> int:
    """Infer completed production steps from JSON/state CSV.

    source:
      - json: use JSON only
      - state: use md_state.csv only
      - max: max(json, state) for robustness to stale single-source metadata
    """
    j_steps = 0
    s_steps = 0

    if json_path is not None and json_path.is_file():
        try:
            j = json.loads(json_path.read_text())
            j_steps = int(j.get("md_production_steps_completed") or j.get("md_production_steps") or 0)
        except Exception:
            j_steps = 0

    if fallback_state_csv is not None and fallback_state_csv.exists():
        try:
            sdf = pd.read_csv(fallback_state_csv)
            for col in ('#"Step"', "Step"):
                if col in sdf.columns:
                    steps = pd.to_numeric(sdf[col], errors="coerce").dropna()
                    if not steps.empty:
                        s_steps = int(steps.max())
                        break
        except Exception:
            s_steps = 0

    src = str(source).strip().lower()
    if src == "json":
        return j_steps
    if src == "state":
        return s_steps
    return max(j_steps, s_steps)


def _load_fold_map() -> dict:
    """Load DOR fold-resistance values from the authoritative xlsx source.

    Falls back to manifests/md_manifest.csv if the xlsx is unavailable.
    """
    xl_path = REPO / "data" / "DRM-susceptibilities.csv.xlsx"
    if xl_path.exists():
        try:
            xl = pd.read_excel(
                xl_path, header=None,
                names=["mutation_raw", "rpv_fold", "dor_fold"],
            )
            # Row 0 is a text header ("Mutations", "RPV...", "DOR...") — drop it.
            xl = xl[xl["mutation_raw"].notna()]
            xl = xl[xl["mutation_raw"].astype(str).str.strip() != "Mutations"]
            fold: dict = {"WT": 1.0}
            for _, row in xl.iterrows():
                mut_raw = str(row["mutation_raw"]).strip()
                # Normalize: "K103N, M230L" → "K103N+M230L"
                mut_norm = re.sub(r",\s*", "+", mut_raw)
                v = row["dor_fold"]
                try:
                    fold[mut_norm] = float(v) if pd.notna(v) else None
                except (TypeError, ValueError):
                    pass
            return fold
        except Exception as exc:
            print(f"  WARNING: could not load xlsx ({exc}), falling back to manifest")
    # Fallback: manifest
    mf = pd.read_csv(REPO / "manifests" / "md_manifest.csv")
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

    out_json_raw = str(row.get("output_json", "")).strip()
    out_json = _remap(Path(out_json_raw)) if out_json_raw else None
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

    rep_dir = topo.parent if topo.exists() else Path(".")
    state_csv = rep_dir / f"{mut_key}_rep{replicate:02d}_md_state.csv"
    total_ns = _infer_total_ns(out_json, fallback_state_csv=state_csv)
    if total_ns is None or not np.isfinite(total_ns) or total_ns <= 0:
        print(f"  WARNING ({mutation} rep{replicate}): could not infer production duration; skipping")
        return []

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
        # Canonical mapping for this project: derive time from production length
        # and frame index; never trust DCD timestamp metadata.
        time_ns = ((float(ts.frame) + 1.0) * float(total_ns)) / max(1, n_frames)
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

    parser = argparse.ArgumentParser(description="Compute DOR torsion-angle traces from MD trajectories.")
    parser.add_argument(
        "--min-production-steps",
        type=int,
        default=0,
        help="Keep only replicates with at least this many completed production steps.",
    )
    parser.add_argument(
        "--step-source",
        choices=["json", "state", "max"],
        default="max",
        help="Metadata source used for min-step filtering (default: max).",
    )
    parser.add_argument(
        "--tag",
        type=str,
        default="",
        help="Optional output tag suffix (e.g. '100ns' -> dor_torsions_100ns.csv / *_100ns.png).",
    )
    args = parser.parse_args()

    run_df = collect_md_results(REPO / "manifests" / "md_manifest.csv")
    FOLD = _load_fold_map()

    if args.min_production_steps > 0:
        keep_idx: list[int] = []
        for i, row in run_df.iterrows():
            mut_key = str(row["safe_label"])
            replicate = int(row["replicate"])
            out_json_raw = str(row.get("output_json", "")).strip()
            out_json = _remap(Path(out_json_raw)) if out_json_raw else None
            topo = _remap(Path(str(row.get("analysis_topology_pdb", "")).strip()))
            rep_dir = topo.parent if topo.exists() else Path(".")
            state_csv = rep_dir / f"{mut_key}_rep{replicate:02d}_md_state.csv"
            done_steps = _infer_completed_steps(
                out_json,
                fallback_state_csv=state_csv,
                source=args.step_source,
            )
            if done_steps >= int(args.min_production_steps):
                keep_idx.append(i)
        run_df = run_df.loc[keep_idx].copy()
        print(
            f"Filtering by completed steps >= {args.min_production_steps} "
            f"(source={args.step_source}): kept {len(run_df)} replicates"
        )

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
    out_csv = OUT_CSV
    tag = args.tag.strip()
    if tag:
        out_csv = REPO / "results" / f"dor_torsions_{tag}.csv"
    df.to_csv(out_csv, index=False)
    print(f"\nWrote {out_csv} ({len(df)} rows)")

    plot(df, FOLD, tag=tag)


# ── plotting ──────────────────────────────────────────────────────────────────
def plot(df: pd.DataFrame, FOLD: dict | None = None, *, tag: str = "") -> None:
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
    suffix = f"_{tag}" if tag else ""
    out = PLOTS_DIR / f"dor_torsions_by_mutation{suffix}.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}")

    # ── Plots 2a-d: one timeseries figure per torsion angle ──
    # X-axis is per-panel and adapts to however much data each mutation has,
    # so partial runs (< 100 ns) are shown as-is without empty trailing space.
    ncols = 5
    nrows = int(np.ceil(n_mut / ncols))
    REP_COLORS = ["#66C2A5", "#FC8D62", "#8DA0CB"]  # ColorBrewer Set2

    TAU_META = [
        ("tau1", "τ₁", "pyridinone–O ether"),
        ("tau2", "τ₂", "O–chlorocyanobenzene ether"),
        ("tau3", "τ₃", "pyridinone-N–CH₂ linker"),
        ("tau4", "τ₄", "CH₂–triazolone linker"),
    ]

    from matplotlib.lines import Line2D
    legend_handles = [
        Line2D([0], [0], color=REP_COLORS[0], lw=2, label="Replicate 1"),
        Line2D([0], [0], color=REP_COLORS[1], lw=2, label="Replicate 2"),
        Line2D([0], [0], color=REP_COLORS[2], lw=2, label="Replicate 3"),
    ]

    for tau_col, tau_label, tau_desc in TAU_META:
        fig, axes = plt.subplots(
            nrows, ncols,
            figsize=(ncols * 3.6, nrows * 3.0),
            squeeze=False,
            sharey=True,
            sharex=False,  # per-panel x so partial runs show their actual range
        )

        for mi, mut in enumerate(mutations):
            row_i = mi // ncols
            col_i = mi % ncols
            ax = axes[row_i, col_i]

            sub = df[df["mutation"] == mut]
            fold = FOLD.get(mut)

            # Per-panel x range based on actual data
            mut_max_t = float(sub["time_ns"].max()) if not sub.empty else 100.0

            # Replicate traces
            for rep_idx, (_, grp) in enumerate(sorted(sub.groupby("replicate"))):
                grp_s = grp.sort_values("time_ns")
                ax.plot(
                    grp_s["time_ns"], grp_s[tau_col],
                    lw=0.9, alpha=0.65,
                    color=REP_COLORS[rep_idx % len(REP_COLORS)],
                    rasterized=True,
                )

            # Reference lines
            ax.axhline(0,   color="#bbbbbb", lw=0.8, zorder=0)
            ax.axhline( 90, color="#e5e5e5", lw=0.5, ls=":", zorder=0)
            ax.axhline(-90, color="#e5e5e5", lw=0.5, ls=":", zorder=0)

            # Title: mutation name + fold change
            if fold is not None and np.isfinite(float(fold)):
                fold_str = "1× (ref)" if mut == "WT" else f"{fold:.0f}× DOR"
            else:
                fold_str = "—"
            ax.set_title(f"{mut}\n{fold_str}", fontsize=8, fontweight="bold",
                         pad=3, linespacing=1.35)

            # Axes decoration
            ax.set_ylim(-185, 185)
            ax.set_xlim(0, mut_max_t)
            ax.set_yticks([-90, 0, 90])
            # Adaptive x ticks: 3 evenly spaced, rounded to nearest 10 ns
            step = max(10, round(mut_max_t / 2 / 10) * 10)
            ax.set_xticks([0, step, min(2 * step, mut_max_t)])
            ax.tick_params(labelsize=7)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            if col_i != 0:
                ax.spines["left"].set_alpha(0.25)
                ax.tick_params(left=False)

            if col_i == 0:
                ax.set_ylabel(f"{tau_label} (°)", fontsize=8)
            if row_i == nrows - 1:
                ax.set_xlabel("Time (ns)", fontsize=8)

        # Hide unused panels
        for mi in range(len(mutations), nrows * ncols):
            axes[mi // ncols][mi % ncols].set_visible(False)

        fig.legend(handles=legend_handles, loc="upper center", ncol=3,
                   fontsize=8.5, frameon=False, bbox_to_anchor=(0.5, 1.02))

        fig.suptitle(
            f"Doravirine {tau_label} torsion ({tau_desc}) — sorted by DOR fold resistance",
            fontsize=10, fontweight="bold", y=1.055,
        )
        fig.tight_layout(h_pad=0.6, w_pad=0.4)

        tau_num = tau_col[-1]  # "1", "2", "3", "4"
        out_tau = PLOTS_DIR / f"dor_torsions_tau{tau_num}_timeseries{suffix}.png"
        fig.savefig(out_tau, dpi=180, bbox_inches="tight")
        plt.close(fig)
        print(f"Wrote {out_tau}")

    # Keep dor_torsions_linker_timeseries.png as a symlink alias for tau4
    import shutil
    tau4_src = PLOTS_DIR / f"dor_torsions_tau4_timeseries{suffix}.png"
    tau4_alias = PLOTS_DIR / f"dor_torsions_linker_timeseries{suffix}.png"
    if tau4_src.exists():
        shutil.copy2(tau4_src, tau4_alias)

    # Remove deprecated legacy panel names to avoid stale, inconsistent outputs.
    legacy_plots = [
        PLOTS_DIR / "dor_torsions_ether_timeseries.png",
    ]
    for legacy in legacy_plots:
        try:
            if legacy.exists():
                legacy.unlink()
        except Exception:
            pass

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
