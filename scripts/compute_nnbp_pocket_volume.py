#!/usr/bin/env python3
"""Compute NNBP pocket volume for both apo and holo simulations.

Pocket center is defined each frame as the centroid of Cα atoms from key
NNBP-lining residues (canonical HIV-1 RT numbering, resid_offset=-3 applied):
  L100, K101, K103, V106, Y181, Y188, F227, P225, W229, L234, H235

This requires no drug to be present and works identically for apo and holo.
Pocket volume = number of grid voxels within POCKET_RADIUS of that centroid
               that are not occluded by any protein heavy atom (VdW excluded).
Uses scipy cKDTree for fast per-frame computation.

Outputs:
  results/nnbp_pocket_volume_dynamics.csv  (per-frame, apo + holo)
  results/plots/nnbp_pocket_volume_apo_vs_holo.png
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
OUT_CSV   = REPO / "results" / "nnbp_pocket_volume_dynamics.csv"
PLOTS_DIR = REPO / "results" / "plots"

# ── pocket geometry ──────────────────────────────────────────────────────────
RESID_OFFSET  = -3
POCKET_RADIUS = 10.0   # Å — generous to capture full NNBP
GRID_SPACING  = 0.75   # Å

# NNBP-lining residues (canonical RT numbering → add RESID_OFFSET for topology)
NNBP_RESIDUES = {
    "L100": 100, "K101": 101, "K103": 103, "V106": 106,
    "Y181": 181, "Y188": 188,
    "F227": 227, "P225": 225, "W229": 229,
    "L234": 234, "H235": 235,
}
# VdW radii for protein heavy atoms
VDW = {"C": 1.7, "N": 1.55, "O": 1.52, "S": 1.8, "P": 1.8}
PROBE_RADIUS  = 1.4   # Å — solvent probe (water)

MUTATION_LABELS = {
    "wt": "WT", "f227c": "F227C", "a98g_f227c": "A98G+F227C",
    "k103n_m230l": "K103N+M230L", "v106i_f227c": "V106I+F227C",
    "v106a": "V106A", "v106a_p225h": "V106A+P225H",
    # holo manifest uses mixed case safe_labels
    "F227C": "F227C", "A98G_F227C": "A98G+F227C",
    "K103N_M230L": "K103N+M230L", "V106I_F227C": "V106I+F227C",
    "V106A": "V106A", "V106A_P225H": "V106A+P225H",
}


# ── helpers ──────────────────────────────────────────────────────────────────
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
        mapped = REPO / s.split(marker, 1)[1]
        if mapped.exists():
            return mapped
    return path


def get_traj_paths(row):
    data = json.loads(Path(row["output_json"]).read_text())
    topo = remap(Path(str(data.get("analysis_topology_pdb", "")).strip()))
    dcd  = remap(Path(str(data.get("analysis_dcd",          "")).strip()))
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


def _build_grid(center: np.ndarray) -> np.ndarray:
    """Return grid points within POCKET_RADIUS of center."""
    s = GRID_SPACING
    r = POCKET_RADIUS
    axes = [np.arange(center[i] - r, center[i] + r + s, s) for i in range(3)]
    grid = np.array(np.meshgrid(*axes, indexing="ij")).reshape(3, -1).T
    grid = grid[np.linalg.norm(grid - center, axis=1) <= r]
    return grid


def pocket_volume_frame(ca_positions: np.ndarray, receptor_pos: np.ndarray,
                        receptor_vdw: np.ndarray) -> float:
    """Compute pocket volume (Å³) for one frame."""
    from scipy.spatial import cKDTree

    center = ca_positions.mean(axis=0)
    grid   = _build_grid(center)
    if len(grid) == 0:
        return 0.0

    tree = cKDTree(receptor_pos)
    # For each grid point find nearest receptor atom; occupied if d < VdW + probe
    # Use a generous fixed exclusion to avoid per-atom radius iteration
    # Then refine with per-atom VdW
    max_vdw = float(receptor_vdw.max()) + PROBE_RADIUS
    candidate_idx = tree.query_ball_point(grid, r=max_vdw)

    free_mask = np.ones(len(grid), dtype=bool)
    for i, neighbours in enumerate(candidate_idx):
        if not neighbours:
            continue
        nb = np.array(neighbours)
        d  = np.linalg.norm(grid[i] - receptor_pos[nb], axis=1)
        if np.any(d < receptor_vdw[nb] + PROBE_RADIUS):
            free_mask[i] = False

    return float(np.sum(free_mask) * GRID_SPACING**3)


def process_replicate(row, leg: str) -> list[dict]:
    import MDAnalysis as mda

    mut_key   = str(row["safe_label"])
    mutation  = MUTATION_LABELS.get(mut_key, mut_key.upper())
    replicate = int(row["replicate"])

    topo, dcd = get_traj_paths(row)
    if not topo.exists() or not dcd.exists():
        print(f"  MISSING: {mutation} rep{replicate} ({leg})")
        return []

    total_ns = infer_total_ns(Path(row["output_json"]))
    u = mda.Universe(str(topo), str(dcd))
    n_frames = len(u.trajectory)

    # Identify p66 (larger subunit) by the segment with the most Cα atoms
    from collections import Counter
    prot_ca_all = u.select_atoms("protein and name CA")
    seg_cnt = Counter(prot_ca_all.segids.tolist())
    p66_seg = max(seg_cnt, key=seg_cnt.get) if seg_cnt else None
    seg_filter = f" and segid {p66_seg}" if p66_seg else ""

    # Select NNBP Cα atoms for pocket center (p66 only)
    resids = [pos + RESID_OFFSET for pos in NNBP_RESIDUES.values()]
    sel_str = ("protein and name CA" + seg_filter +
               " and (" + " or ".join(f"resid {r}" for r in resids) + ")")
    ca_sel = u.select_atoms(sel_str)
    if ca_sel.n_atoms < 5:
        print(f"  WARNING: only {ca_sel.n_atoms} NNBP Cα found for {mutation} rep{replicate} (seg={p66_seg})")

    # Pre-build receptor (all protein heavy atoms, no H)
    receptor = u.select_atoms("protein and not name H*")

    # Pre-compute VdW radii for receptor atoms (static per topology)
    rec_vdw = np.array([
        VDW.get((a.element or a.name[0]).upper()[:1], 1.7)
        for a in receptor.atoms
    ], dtype=float)

    rows_out = []
    for ts in u.trajectory:
        time_ns = (float(ts.frame) / max(1, n_frames - 1)) * total_ns
        vol = pocket_volume_frame(ca_sel.positions.copy(),
                                  receptor.positions.copy(),
                                  rec_vdw)
        rows_out.append({
            "leg":       leg,
            "mutation":  mutation,
            "safe_label": mut_key,
            "replicate": replicate,
            "frame":     int(ts.frame),
            "time_ns":   time_ns,
            "pocket_volume_A3": vol,
        })

    vols = [r["pocket_volume_A3"] for r in rows_out]
    print(f"  {leg:4s} {mutation} rep{replicate}: "
          f"mean={np.mean(vols):.0f}  std={np.std(vols):.0f}  "
          f"min={np.min(vols):.0f}  max={np.max(vols):.0f}  Å³")
    return rows_out


# ── main ─────────────────────────────────────────────────────────────────────
def main():
    # Mutations that have both apo and holo
    apo_mf  = pd.read_csv(REPO / "results" / "apo_md_manifest.csv")
    holo_mf = pd.read_csv(REPO / "results" / "md_manifest.csv")

    apo_labels  = set(apo_mf["safe_label"].str.lower().unique())
    holo_subset = holo_mf[holo_mf["safe_label"].str.lower().isin(apo_labels)].copy()

    print(f"Apo  replicates: {len(apo_mf)}")
    print(f"Holo replicates (matched): {len(holo_subset)}")
    print(f"Mutations: {sorted(apo_labels)}\n")

    all_rows: list[dict] = []

    print("── APO ──")
    for _, row in apo_mf.iterrows():
        all_rows.extend(process_replicate(row, "apo"))

    print("\n── HOLO ──")
    for _, row in holo_subset.iterrows():
        all_rows.extend(process_replicate(row, "holo"))

    if not all_rows:
        print("No data collected.")
        return

    df = pd.DataFrame(all_rows)
    df.to_csv(OUT_CSV, index=False)
    print(f"\nWrote {OUT_CSV} ({len(df)} rows)")

    FOLD = _load_fold_map()
    plot(df, FOLD)


# ── plotting ──────────────────────────────────────────────────────────────────
def plot(df: pd.DataFrame, FOLD: dict | None = None):
    if FOLD is None:
        FOLD = _load_fold_map()
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    mutations_ordered = sorted(
        df["mutation"].unique(),
        key=lambda m: FOLD.get(m, 999)
    )
    n = len(mutations_ordered)
    fig, axes = plt.subplots(1, n, figsize=(3.2 * n, 4.5), sharey=True)
    if n == 1:
        axes = [axes]

    for ax, mut in zip(axes, mutations_ordered):
        apo_vals  = df[(df["mutation"] == mut) & (df["leg"] == "apo" )]["pocket_volume_A3"].to_numpy()
        holo_vals = df[(df["mutation"] == mut) & (df["leg"] == "holo")]["pocket_volume_A3"].to_numpy()

        if len(apo_vals) == 0 or len(holo_vals) == 0:
            ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)
            continue
        vp = ax.violinplot([apo_vals, holo_vals], positions=[0, 1],
                           showmedians=True, showextrema=False)
        vp["bodies"][0].set_facecolor("#5B9BD5"); vp["bodies"][0].set_alpha(0.7)
        vp["bodies"][1].set_facecolor("#ED7D31"); vp["bodies"][1].set_alpha(0.7)
        vp["cmedians"].set_color("black"); vp["cmedians"].set_linewidth(1.5)

        ax.set_xticks([0, 1])
        ax.set_xticklabels(["apo", "holo"], fontsize=8)
        fold = FOLD.get(mut, "?")
        ax.set_title(f"{mut}\n({fold:.0f}×)", fontsize=8, fontweight="bold")
        ax.grid(axis="y", alpha=0.25, linestyle=":")

        # Annotate medians
        for x, vals, color in [(0, apo_vals, "#1F497D"), (1, holo_vals, "#843C0C")]:
            if len(vals):
                ax.text(x, np.median(vals) + 15, f"{np.median(vals):.0f}",
                        ha="center", va="bottom", fontsize=7, color=color, fontweight="bold")

    axes[0].set_ylabel("NNBP pocket volume (Å³)")
    fig.suptitle("NNBP pocket volume: apo vs holo\n"
                 "(center = NNBP residue Cα centroid, radius=10 Å, probe=1.4 Å)",
                 fontsize=10, fontweight="bold")
    # Legend
    from matplotlib.patches import Patch
    axes[-1].legend(handles=[Patch(facecolor="#5B9BD5", alpha=0.7, label="apo"),
                              Patch(facecolor="#ED7D31", alpha=0.7, label="holo")],
                    fontsize=8, frameon=False, loc="upper right")
    fig.tight_layout()
    out = PLOTS_DIR / "nnbp_pocket_volume_apo_vs_holo.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}")

    # Summary table
    print("\n─── Pocket volume summary (median ± IQR/2, Å³) ───")
    print(f"{'Mutation':<20} {'Fold':>6}  {'apo median':>12}  {'holo median':>12}  {'Δ (holo-apo)':>14}")
    for mut in mutations_ordered:
        a = df[(df["mutation"] == mut) & (df["leg"] == "apo" )]["pocket_volume_A3"].to_numpy()
        h = df[(df["mutation"] == mut) & (df["leg"] == "holo")]["pocket_volume_A3"].to_numpy()
        if len(a) and len(h):
            delta = np.median(h) - np.median(a)
            print(f"{mut:<20} {FOLD.get(mut,'?'):>6.0f}  "
                  f"{np.median(a):>10.0f} Å³  {np.median(h):>10.0f} Å³  "
                  f"{delta:>+12.0f} Å³")


if __name__ == "__main__":
    main()
