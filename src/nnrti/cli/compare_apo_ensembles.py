#!/usr/bin/env python3
"""WT-referenced comparison of apo pocket ensembles, on a matched time window.

Reviewer 1 (comment 2) asked for longer apo simulations of selected mutants,
P225- and M230-containing in particular, compared against WT, on the grounds
that the NNIBP is induced and mutations may alter the apo conformation before
DOR ever binds.

Apo coverage on disk is uneven: three replicates at 100 ns for A98G+F227C,
F227C and G190E, 37-75 ns for K103N+M230L, and 10 ns for WT, V106A,
V106A+P225H and V106I+F227C. WT is the reference, so its 10 ns caps any
like-for-like comparison. This analysis therefore truncates every trajectory to
its first 10 ns and compares on that matched window, which is what the existing
data supports without new simulation.

It also asks the separate question that decides whether new simulation is
needed at all: for the genotypes that do run to 100 ns, does the observable
drift after the first 10 ns? If the long runs are flat beyond 10 ns, the short
WT reference is a minor limitation and can be stated as one. If they drift, the
drift tells us how long a WT extension would have to be.

Observables are those that would report pocket collapse: the NNIBP volume, the
chi1 rotamers of Tyr181 and Tyr188 (whose side chains occupy the pocket in
unliganded RT), and the Ca RMSF of the pocket-lining residues.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.compare_apo_ensembles
"""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

RESID_OFFSET = -3
NNIBP_P66 = (100, 101, 103, 106, 107, 108, 179, 181, 188, 189, 190, 227, 229, 234, 318)
MATCH_NS = 10.0


def chi1(u, resid: int, segid: str) -> float:
    """Chi1 dihedral (N-CA-CB-CG/OG/SG) in degrees, or nan."""
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import calc_dihedrals
    sel = u.select_atoms(f"segid {segid} and resid {resid}")
    names = set(sel.names)
    g = next((n for n in ("CG", "CG1", "OG", "OG1", "SG", "CG2") if n in names), None)
    if g is None or not {"N", "CA", "CB"} <= names:
        return float("nan")
    p = [sel.select_atoms(f"name {n}").positions[0] for n in ("N", "CA", "CB", g)]
    return float(np.degrees(calc_dihedrals(p[0][None], p[1][None], p[2][None], p[3][None])[0]))


def main(argv: list[str] | None = None) -> int:
    import MDAnalysis as mda
    import mdtraj as md
    from nnrti.analysis.metrics import pocket_volume_proxy_from_universe, _p66_segid

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--apo-root", type=Path, default=Path("results/md_runs/apo"))
    p.add_argument("--output-dir", type=Path, default=Path("results/analysis/apo_ensembles"))
    p.add_argument("--match-ns", type=float, default=MATCH_NS)
    args = p.parse_args(argv)

    runs = []
    for d in sorted(args.apo_root.glob("*/rep_*")):
        dcd = sorted(d.glob("*analysis.dcd"))
        top = sorted(d.glob("*analysis_topology.pdb"))
        js = sorted(d.glob("*_apo_rep*.json"))
        if not (dcd and top and js):
            continue
        import json
        meta = json.loads(js[0].read_text())
        steps = meta.get("md_production_steps_completed")
        if not steps:
            continue
        runs.append((d.parent.name, int(d.name.split("_")[-1]), dcd[0], top[0],
                     steps * 2e-6))          # 2 fs timestep -> ns
    if not runs:
        raise SystemExit("no usable apo runs found")

    rows = []
    for label, rep, dcd, top, total_ns in runs:
        t = md.load(str(dcd), top=str(top))
        dt = total_ns / t.n_frames
        n_match = max(1, int(round(args.match_ns / dt)))
        u = mda.Universe(str(top))
        seg = _p66_segid(u)
        ca = u.select_atoms(f"segid {seg} and name CA and resid "
                            + " ".join(str(r + RESID_OFFSET) for r in NNIBP_P66))
        traj = t.xyz * 10.0                   # nm -> A

        for i in range(t.n_frames):
            u.atoms.positions = traj[i]
            window = "matched_10ns" if i < n_match else "beyond_10ns"
            rows.append({
                "safe_label": label, "replicate": rep, "frame": i,
                "time_ns": (i + 1) * dt, "window": window,
                "total_ns": total_ns,
                "pocket_volume": pocket_volume_proxy_from_universe(u, ""),
                "chi1_Y181": chi1(u, 181 + RESID_OFFSET, seg),
                "chi1_Y188": chi1(u, 188 + RESID_OFFSET, seg),
            })
        # pocket Ca RMSF over the matched window
        sub = t[:n_match]
        sub.superpose(sub, 0)
        idx = [a.index for a in sub.topology.atoms
               if a.name == "CA" and (a.residue.resSeq in [r + RESID_OFFSET for r in NNIBP_P66])]
        if idx:
            xyz = sub.xyz[:, idx, :] * 10.0
            rmsf = float(np.sqrt(((xyz - xyz.mean(0)) ** 2).sum(-1).mean(0)).mean())
        else:
            rmsf = float("nan")
        rows[-1]["pocket_ca_rmsf_matched"] = rmsf
        print(f"  {label:14s} rep{rep}  {t.n_frames:4d} frames  {total_ns:6.1f} ns  "
              f"dt={dt:.3f} ns  matched={n_match} frames  RMSF={rmsf:.2f} A", flush=True)

    d = pd.DataFrame(rows)
    out = args.output_dir / "tables"
    out.mkdir(parents=True, exist_ok=True)
    d.to_csv(out / "apo_frame_metrics.csv", index=False)

    m = d[d.window == "matched_10ns"]
    summ = (m.groupby("safe_label")
              .agg(n_rep=("replicate", "nunique"),
                   pocket_volume=("pocket_volume", "mean"),
                   pocket_volume_sd=("pocket_volume", "std"),
                   chi1_Y181=("chi1_Y181", "mean"),
                   chi1_Y188=("chi1_Y188", "mean"))
              .round(1))
    if "wt" in summ.index:
        summ["dVolume_vs_WT"] = (summ.pocket_volume - summ.loc["wt", "pocket_volume"]).round(1)
    summ.to_csv(out / "apo_matched_window_summary.csv")

    print("\n=== matched 10 ns window, apo ===")
    print(summ.to_string())

    drift = []
    for lbl, g in d[d.total_ns > 20].groupby("safe_label"):
        a = g[g.window == "matched_10ns"].pocket_volume.mean()
        b = g[g.window == "beyond_10ns"].pocket_volume.mean()
        drift.append({"safe_label": lbl, "first_10ns": round(a, 1),
                      "beyond_10ns": round(b, 1), "drift": round(b - a, 1)})
    if drift:
        dr = pd.DataFrame(drift)
        dr.to_csv(out / "apo_drift_beyond_10ns.csv", index=False)
        print("\n=== does the pocket volume drift after 10 ns? (long runs only) ===")
        print(dr.to_string(index=False))
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
