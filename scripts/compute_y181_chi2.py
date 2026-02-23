#!/usr/bin/env python3
"""Compute Y181 chi2 dihedral from apo MD trajectories.

Y181 chi2 = CA-CB-CG-CD1 dihedral.

Crystal structure reference values (p66 chain):
  4NCG (DOR-bound, NNBP open)  : chi2 = +79°
  1DLO (apo, open-like)        : chi2 = +87°
  1HMV (apo, NNBP closed)      : chi2 = -18°

The ~95° difference in chi2 distinguishes open (ring pointing out into NNBP)
from closed (ring rotated inward, pocket collapsed).
"""
import json
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[1]
MANIFEST = REPO / "results" / "apo_md_manifest.csv"
OUT_CSV = REPO / "results" / "apo_y181_chi2.csv"
PLOTS_DIR = REPO / "results" / "plots"

RESID_OFFSET = -3
Y181_RESID = 181 + RESID_OFFSET  # = 178

# Reference chi2 values from crystal structures
CHI2_OPEN   =  79.0   # 4NCG (DOR-bound p66)
CHI2_CLOSED = -18.0   # 1HMV (apo closed p66)
# Threshold: < +30° considered "closed-like"
CHI2_CLOSED_THRESH = 30.0

MUTATION_LABELS = {
    "wt": "WT", "f227c": "F227C", "a98g_f227c": "A98G+F227C",
    "k103n_m230l": "K103N+M230L", "v106i_f227c": "V106I+F227C",
    "v106a": "V106A", "v106a_p225h": "V106A+P225H",
}


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


def remap(path: Path) -> Path:
    if path.exists():
        return path
    marker = "nnrti-mechanisms/"
    s = str(path)
    if marker in s:
        rel = s.split(marker, 1)[1]
        mapped = REPO / rel
        if mapped.exists():
            return mapped
    return path


def get_traj_paths(row):
    data = json.loads(Path(row["output_json"]).read_text())
    topo = remap(Path(str(data.get("analysis_topology_pdb", "")).strip()))
    dcd  = remap(Path(str(data.get("analysis_dcd", "")).strip()))
    return topo, dcd


def infer_total_ns(json_path: Path) -> float:
    try:
        j = json.loads(json_path.read_text())
        steps = int(j.get("md_production_steps_completed") or j.get("md_production_steps") or 0)
        if steps > 0:
            return steps * 2.0 / 1_000_000.0
    except Exception:
        pass
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
    return 10.0


def calc_dihedral(p0, p1, p2, p3):
    b0 = p0 - p1; b1 = p2 - p1; b2 = p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    return np.degrees(np.arctan2(np.dot(np.cross(b1n, v), w), np.dot(v, w)))


def process_replicate(row):
    import MDAnalysis as mda

    mut_key  = str(row["safe_label"])
    mutation = MUTATION_LABELS.get(mut_key, mut_key.upper())
    replicate = int(row["replicate"])

    topo, dcd = get_traj_paths(row)
    if not topo.exists() or not dcd.exists():
        print(f"  MISSING: {mutation} rep{replicate}")
        return []

    total_ns = infer_total_ns(Path(row["output_json"]))
    u = mda.Universe(str(topo), str(dcd))
    n_frames = len(u.trajectory)

    sel_ca  = u.select_atoms(f"protein and resid {Y181_RESID} and name CA")
    sel_cb  = u.select_atoms(f"protein and resid {Y181_RESID} and name CB")
    sel_cg  = u.select_atoms(f"protein and resid {Y181_RESID} and name CG")
    sel_cd1 = u.select_atoms(f"protein and resid {Y181_RESID} and name CD1")

    missing = [n for n, ag in [("CA", sel_ca), ("CB", sel_cb),
                                ("CG", sel_cg), ("CD1", sel_cd1)] if ag.n_atoms == 0]
    if missing:
        print(f"  MISSING Y181 atoms {missing} for {mutation} rep{replicate} (resid {Y181_RESID})")
        for r in range(Y181_RESID - 2, Y181_RESID + 3):
            test = u.select_atoms(f"protein and resid {r} and resname TYR")
            if test.n_atoms > 0:
                print(f"    TYR found at resid {r}")
        return []

    rows = []
    for ts in u.trajectory:
        time_ns = (float(ts.frame) / max(1, n_frames - 1)) * total_ns
        chi2 = calc_dihedral(
            sel_ca.positions[0], sel_cb.positions[0],
            sel_cg.positions[0], sel_cd1.positions[0],
        )
        rows.append({
            "mutation": mutation, "safe_label": mut_key,
            "replicate": replicate, "frame": int(ts.frame),
            "time_ns": time_ns, "y181_chi2_deg": chi2,
        })

    vals = [r["y181_chi2_deg"] for r in rows]
    print(f"  {mutation} rep{replicate}: {len(rows)} frames, "
          f"chi2 mean={np.mean(vals):.1f}°  std={np.std(vals):.1f}°  "
          f"frac_closed={np.mean(np.array(vals) < CHI2_CLOSED_THRESH):.3f}")
    return rows


def main():
    mf = pd.read_csv(MANIFEST)
    print(f"Processing {len(mf)} apo replicates for Y181 chi2...")
    print(f"Reference: open={CHI2_OPEN}° (4NCG), closed={CHI2_CLOSED}° (1HMV), threshold={CHI2_CLOSED_THRESH}°\n")

    all_rows = []
    for _, row in mf.iterrows():
        all_rows.extend(process_replicate(row))

    if not all_rows:
        print("No data collected.")
        return

    df = pd.DataFrame(all_rows)
    df.to_csv(OUT_CSV, index=False)
    print(f"\nWrote {OUT_CSV} ({len(df)} rows)")

    FOLD = _load_fold_map()
    plot(df, FOLD)


def plot(df, FOLD=None):
    if FOLD is None:
        FOLD = _load_fold_map()
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    mutations_ordered = sorted(df["mutation"].unique(), key=lambda m: FOLD.get(m, 999))
    n_mut = len(mutations_ordered)
    colors = dict(zip(mutations_ordered, cm.RdYlGn_r(np.linspace(0.1, 0.9, n_mut))))

    # --- Timeseries ---
    fig, axes = plt.subplots(n_mut, 1, figsize=(12, 2.2 * n_mut), sharex=False)
    if n_mut == 1:
        axes = [axes]

    for ax, mut in zip(axes, mutations_ordered):
        sub = df[df["mutation"] == mut]
        for rep, grp in sub.groupby("replicate"):
            g = grp.sort_values("time_ns")
            ax.plot(g["time_ns"], g["y181_chi2_deg"], alpha=0.7, linewidth=0.8,
                    color=colors[mut], label=f"rep{rep}")
        ax.axhline(CHI2_OPEN,   color="darkgreen", linestyle="--", linewidth=1.0, alpha=0.8,
                   label=f"open 4NCG ({CHI2_OPEN:.0f}°)")
        ax.axhline(CHI2_CLOSED, color="navy",      linestyle="--", linewidth=1.0, alpha=0.8,
                   label=f"closed 1HMV ({CHI2_CLOSED:.0f}°)")
        ax.axhline(CHI2_CLOSED_THRESH, color="gray", linestyle=":", linewidth=0.8, alpha=0.6,
                   label=f"threshold ({CHI2_CLOSED_THRESH:.0f}°)")
        ax.set_ylabel("χ2 (°)", fontsize=8)
        ax.set_ylim(-185, 185)
        ax.set_yticks([-180, -90, -18, 0, 30, 79, 90, 180])
        fold = FOLD.get(mut, "?")
        ax.set_title(f"{mut}  (fold={fold:.0f}×)", fontsize=9, loc="left", fontweight="bold")
        ax.grid(alpha=0.2, linestyle=":")
        ax.legend(fontsize=6, ncol=4, frameon=False, loc="upper right")

    axes[-1].set_xlabel("Time (ns)")
    fig.suptitle("Y181 χ2 dihedral — apo simulations\n"
                 f"(open={CHI2_OPEN:.0f}° [4NCG], closed={CHI2_CLOSED:.0f}° [1HMV])",
                 fontsize=11, fontweight="bold")
    fig.tight_layout()
    out = PLOTS_DIR / "apo_y181_chi2_timeseries.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}")

    # --- Distribution ---
    fig, ax = plt.subplots(figsize=(10, 4))
    bins = np.linspace(-180, 180, 73)
    for mut in mutations_ordered:
        vals = df[df["mutation"] == mut]["y181_chi2_deg"].to_numpy()
        counts, edges = np.histogram(vals, bins=bins, density=True)
        centers = (edges[:-1] + edges[1:]) / 2
        fold = FOLD.get(mut, "?")
        ax.plot(centers, counts, color=colors[mut], linewidth=1.5,
                label=f"{mut} ({fold:.0f}×)")
        ax.fill_between(centers, counts, alpha=0.08, color=colors[mut])

    ax.axvline(CHI2_OPEN,   color="darkgreen", linestyle="--", linewidth=1.2, alpha=0.8,
               label=f"open 4NCG ({CHI2_OPEN:.0f}°)")
    ax.axvline(CHI2_CLOSED, color="navy",      linestyle="--", linewidth=1.2, alpha=0.8,
               label=f"closed 1HMV ({CHI2_CLOSED:.0f}°)")
    ax.axvline(CHI2_CLOSED_THRESH, color="gray", linestyle=":", linewidth=0.8, alpha=0.6,
               label=f"threshold ({CHI2_CLOSED_THRESH:.0f}°)")
    ax.set_xlabel("Y181 χ2 (°)")
    ax.set_ylabel("Density")
    ax.set_title("Y181 χ2 rotamer distribution — apo simulations", fontweight="bold")
    ax.set_xlim(-180, 180)
    ax.set_xticks([-180, -90, -18, 0, 30, 60, 79, 90, 180])
    ax.legend(fontsize=8, frameon=False)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    out = PLOTS_DIR / "apo_y181_chi2_distribution.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}")

    # --- Fraction closed vs fold ---
    records = []
    for mut in mutations_ordered:
        sub = df[df["mutation"] == mut]
        for rep, grp in sub.groupby("replicate"):
            vals = grp["y181_chi2_deg"].to_numpy()
            records.append({
                "mutation": mut, "replicate": rep,
                "frac_closed": float(np.mean(vals < CHI2_CLOSED_THRESH)),
                "fold": FOLD.get(mut, np.nan),
            })
    rdf = pd.DataFrame(records)

    fig, ax = plt.subplots(figsize=(7, 4))
    for mut in mutations_ordered:
        sub = rdf[rdf["mutation"] == mut]
        mean_fc = sub["frac_closed"].mean()
        ax.scatter(sub["fold"], sub["frac_closed"],
                   color=colors[mut], s=50, zorder=3, alpha=0.7)
        ax.scatter(FOLD.get(mut, np.nan), mean_fc,
                   color=colors[mut], s=120, marker="D",
                   edgecolors="black", linewidths=0.8, zorder=4,
                   label=f"{mut} ({FOLD.get(mut,'?'):.0f}×)")
    ax.set_xscale("log")
    ax.set_xlabel("Fold resistance (log scale)")
    ax.set_ylabel(f"Fraction of frames with χ2 < {CHI2_CLOSED_THRESH:.0f}° (closed-like)")
    ax.set_title("Y181 ring closure tendency vs. DOR resistance — apo sims", fontweight="bold")
    ax.legend(fontsize=8, frameon=False)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    out = PLOTS_DIR / "apo_y181_chi2_vs_fold.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}")

    # Summary table
    print(f"\n--- Y181 chi2 summary (frac time chi2 < {CHI2_CLOSED_THRESH}°, closed-like) ---")
    print(f"{'Mutation':<20} {'Fold':>6}  {'n_reps':>6}  {'mean chi2':>10}  {'frac_closed (mean±std)':>25}")
    for mut in mutations_ordered:
        sub_df = df[df["mutation"] == mut]
        sub_r  = rdf[rdf["mutation"] == mut]
        mean_chi2 = sub_df["y181_chi2_deg"].mean()
        m, s = sub_r["frac_closed"].mean(), sub_r["frac_closed"].std()
        print(f"{mut:<20} {FOLD.get(mut,'?'):>6.0f}  {len(sub_r):>6}  "
              f"{mean_chi2:>10.1f}°  {m:.3f} ± {s:.3f}")


if __name__ == "__main__":
    main()
