#!/usr/bin/env python3
"""Compute Y181 chi1 dihedral from apo MD trajectories.

Y181 chi1 = N-CA-CB-CG dihedral of residue Y181 (resid 178 with offset -3).
- Open state (NNBP formed): chi1 ~ +60° (gauche+) or 180° (trans)
- Closed state (Y181 rotated inward): chi1 ~ -60° (gauche-)

Plots timeseries for all 7 apo mutations and a combined distribution.
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
OUT_CSV = REPO / "results" / "apo_y181_chi1.csv"
PLOTS_DIR = REPO / "results" / "plots"

RESID_OFFSET = -3
Y181_RESID = 181 + RESID_OFFSET  # = 178

# Fold resistance for ordering/coloring
FOLD = {
    "WT": 1, "F227C": 2, "A98G+F227C": 4,
    "K103N+M230L": 8, "V106I+F227C": 38,
    "V106A": 41, "V106A+P225H": 153,
}

MUTATION_LABELS = {
    "wt": "WT", "f227c": "F227C", "a98g_f227c": "A98G+F227C",
    "k103n_m230l": "K103N+M230L", "v106i_f227c": "V106I+F227C",
    "v106a": "V106A", "v106a_p225h": "V106A+P225H",
}


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
    dcd = remap(Path(str(data.get("analysis_dcd", "")).strip()))
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


def compute_dihedral(p0, p1, p2, p3):
    """Praxitelous dihedral angle in degrees."""
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2
    b1_norm = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1_norm) * b1_norm
    w = b2 - np.dot(b2, b1_norm) * b1_norm
    x = np.dot(v, w)
    y = np.dot(np.cross(b1_norm, v), w)
    return np.degrees(np.arctan2(y, x))


def process_replicate(row):
    import MDAnalysis as mda

    mut_key = str(row["safe_label"])
    mutation = MUTATION_LABELS.get(mut_key, mut_key.upper())
    replicate = int(row["replicate"])

    topo, dcd = get_traj_paths(row)
    if not topo.exists() or not dcd.exists():
        print(f"  MISSING: {mutation} rep{replicate}")
        return []

    total_ns = infer_total_ns(Path(row["output_json"]))
    u = mda.Universe(str(topo), str(dcd))
    n_frames = len(u.trajectory)

    # Select chi1 atoms: N, CA, CB, CG of Y181 (resid 178)
    sel_n  = u.select_atoms(f"protein and resid {Y181_RESID} and name N")
    sel_ca = u.select_atoms(f"protein and resid {Y181_RESID} and name CA")
    sel_cb = u.select_atoms(f"protein and resid {Y181_RESID} and name CB")
    sel_cg = u.select_atoms(f"protein and resid {Y181_RESID} and name CG")

    if any(ag.n_atoms == 0 for ag in [sel_n, sel_ca, sel_cb, sel_cg]):
        print(f"  MISSING Y181 atoms for {mutation} rep{replicate} (resid {Y181_RESID})")
        # Try nearby resids
        for r in range(Y181_RESID - 2, Y181_RESID + 3):
            test = u.select_atoms(f"protein and resid {r} and resname TYR")
            if test.n_atoms > 0:
                print(f"    TYR found at resid {r}")
        return []

    rows = []
    for ts in u.trajectory:
        time_ns = (float(ts.frame) / max(1, n_frames - 1)) * total_ns
        chi1 = compute_dihedral(
            sel_n.positions[0],
            sel_ca.positions[0],
            sel_cb.positions[0],
            sel_cg.positions[0],
        )
        rows.append({
            "mutation": mutation,
            "safe_label": mut_key,
            "replicate": replicate,
            "frame": int(ts.frame),
            "time_ns": time_ns,
            "y181_chi1_deg": chi1,
        })

    print(f"  {mutation} rep{replicate}: {len(rows)} frames, "
          f"chi1 mean={np.mean([r['y181_chi1_deg'] for r in rows]):.1f}° "
          f"std={np.std([r['y181_chi1_deg'] for r in rows]):.1f}°")
    return rows


def main():
    import MDAnalysis as mda  # noqa: ensure import

    mf = pd.read_csv(MANIFEST)
    print(f"Processing {len(mf)} apo replicates...")

    all_rows = []
    for _, row in mf.iterrows():
        all_rows.extend(process_replicate(row))

    if not all_rows:
        print("No data collected.")
        return

    df = pd.DataFrame(all_rows)
    df.to_csv(OUT_CSV, index=False)
    print(f"\nWrote {OUT_CSV} ({len(df)} rows)")

    plot(df)


def plot(df):
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    mutations_ordered = sorted(
        df["mutation"].unique(),
        key=lambda m: FOLD.get(m, 999)
    )
    n_mut = len(mutations_ordered)
    colors = dict(zip(mutations_ordered, cm.RdYlGn_r(np.linspace(0.1, 0.9, n_mut))))

    # --- Panel 1: chi1 timeseries, one row per mutation ---
    fig, axes = plt.subplots(n_mut, 1, figsize=(12, 2.2 * n_mut), sharex=False)
    if n_mut == 1:
        axes = [axes]

    for ax, mut in zip(axes, mutations_ordered):
        sub = df[df["mutation"] == mut]
        for rep, grp in sub.groupby("replicate"):
            g = grp.sort_values("time_ns")
            ax.plot(g["time_ns"], g["y181_chi1_deg"], alpha=0.7, linewidth=0.8,
                    color=colors[mut], label=f"rep{rep}")
        ax.axhline(-60, color="navy", linestyle="--", linewidth=0.8, alpha=0.6, label="gauche- (closed)")
        ax.axhline(180, color="darkgreen", linestyle="--", linewidth=0.8, alpha=0.6, label="trans (open)")
        ax.axhline(60, color="gray", linestyle=":", linewidth=0.8, alpha=0.6, label="gauche+")
        ax.set_ylabel("χ1 (°)", fontsize=8)
        ax.set_ylim(-185, 185)
        ax.set_yticks([-180, -60, 0, 60, 180])
        fold = FOLD.get(mut, "?")
        ax.set_title(f"{mut}  (fold={fold}×)", fontsize=9, loc="left", fontweight="bold")
        ax.grid(alpha=0.2, linestyle=":")
        ax.legend(fontsize=6, ncol=4, frameon=False, loc="upper right")

    axes[-1].set_xlabel("Time (ns)")
    fig.suptitle("Y181 χ1 dihedral — apo simulations\n"
                 "(gauche− = closed/inward, trans = open/NNBP)", fontsize=11, fontweight="bold")
    fig.tight_layout()
    out = PLOTS_DIR / "apo_y181_chi1_timeseries.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}")

    # --- Panel 2: chi1 distribution (KDE/histogram) per mutation ---
    fig, ax = plt.subplots(figsize=(10, 4))
    bins = np.linspace(-180, 180, 73)
    for mut in mutations_ordered:
        vals = df[df["mutation"] == mut]["y181_chi1_deg"].to_numpy()
        counts, edges = np.histogram(vals, bins=bins, density=True)
        centers = (edges[:-1] + edges[1:]) / 2
        ax.plot(centers, counts, color=colors[mut], linewidth=1.5,
                label=f"{mut} ({FOLD.get(mut,'?')}×)")
        ax.fill_between(centers, counts, alpha=0.08, color=colors[mut])

    ax.axvline(-60, color="navy", linestyle="--", linewidth=1, alpha=0.7, label="gauche- (closed)")
    ax.axvline(180, color="darkgreen", linestyle="--", linewidth=1, alpha=0.7, label="trans (open)")
    ax.axvline(60, color="gray", linestyle=":", linewidth=1, alpha=0.5, label="gauche+")
    ax.set_xlabel("Y181 χ1 (°)")
    ax.set_ylabel("Density")
    ax.set_title("Y181 χ1 rotamer distribution — apo simulations", fontweight="bold")
    ax.set_xlim(-180, 180)
    ax.set_xticks([-180, -120, -60, 0, 60, 120, 180])
    ax.legend(fontsize=8, frameon=False)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    out = PLOTS_DIR / "apo_y181_chi1_distribution.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}")

    # --- Panel 3: fraction of time in closed rotamer (chi1 < -30°) ---
    records = []
    for mut in mutations_ordered:
        sub = df[df["mutation"] == mut]
        for rep, grp in sub.groupby("replicate"):
            vals = grp["y181_chi1_deg"].to_numpy()
            frac_closed = np.mean(vals < -30)
            records.append({"mutation": mut, "replicate": rep,
                            "frac_closed": frac_closed,
                            "fold": FOLD.get(mut, np.nan)})
    rdf = pd.DataFrame(records)

    fig, ax = plt.subplots(figsize=(7, 4))
    for i, mut in enumerate(mutations_ordered):
        sub = rdf[rdf["mutation"] == mut]
        mean_fc = sub["frac_closed"].mean()
        ax.scatter(sub["fold"], sub["frac_closed"], color=colors[mut],
                   s=60, zorder=3, alpha=0.8)
        ax.scatter(FOLD.get(mut, np.nan), mean_fc, color=colors[mut],
                   s=120, marker="D", edgecolors="black", linewidths=0.8, zorder=4,
                   label=mut)

    ax.set_xscale("log")
    ax.set_xlabel("Fold resistance (log scale)")
    ax.set_ylabel("Fraction of frames with χ1 < −30° (closed rotamer)")
    ax.set_title("Y181 closure tendency vs. DOR resistance — apo sims", fontweight="bold")
    ax.legend(fontsize=8, frameon=False)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    out = PLOTS_DIR / "apo_y181_chi1_vs_fold.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}")

    # Summary stats
    print("\n--- Y181 chi1 summary (fraction of time in closed rotamer chi1 < -30°) ---")
    print(f"{'Mutation':<20} {'Fold':>6}  {'n_reps':>6}  {'frac_closed (mean±std)':>25}")
    for mut in mutations_ordered:
        sub = rdf[rdf["mutation"] == mut]
        m, s = sub["frac_closed"].mean(), sub["frac_closed"].std()
        print(f"{mut:<20} {FOLD.get(mut,'?'):>6}  {len(sub):>6}  {m:.3f} ± {s:.3f}")


if __name__ == "__main__":
    main()
