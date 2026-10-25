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

# A replicate that was extended keeps its original frames and appends new ones
# at the longer run's stride, so the analysis DCD is not uniformly spaced.
# Both segments follow interval = production_steps // 200: the original 10 ns
# run gives 200 frames of 0.05 ns, the extension to 100 ns gives 0.5 ns each.
FIRST_SEGMENT_FRAMES = 200
FIRST_SEGMENT_DT_NS = 0.05
EXTENSION_DT_NS = 0.5

# Pocket-lining Ca radius of gyration. Healthy replicates sit at 7.2-7.8 A;
# a replicate whose pocket has come apart reaches 20 A and reports a pocket
# volume in the thousands, which destroys any genotype mean it enters.
POCKET_RG_MAX_A = 10.0


def true_steps_completed(rep_dir: Path, json_steps: int) -> int:
    """Steps actually on disk, from the live state CSV.

    The run JSON records md_production_steps_completed only when a run
    finishes, so an in-flight replicate still advertises its previous target
    -- WT rep_02 and rep_03 claimed 10 ns while holding 91 and 83 ns of
    trajectory. The state CSV is written continuously and does not lie.
    """
    best = int(json_steps or 0)
    for csv_path in sorted(rep_dir.glob("*_md_state.csv")):
        try:
            with csv_path.open() as handle:
                handle.readline()
                for line in handle:
                    line = line.strip()
                    if not line:
                        continue
                    try:
                        best = max(best, int(float(line.split(",")[0].strip('"'))))
                    except ValueError:
                        continue
        except OSError:
            continue
    return best


def frame_times_ns(n_frames: int, total_ns: float) -> np.ndarray:
    """Per-frame simulation time, honouring the two-stride schedule."""
    if n_frames <= FIRST_SEGMENT_FRAMES:
        dt = total_ns / max(1, n_frames)
        return np.arange(1, n_frames + 1) * dt
    head = np.arange(1, FIRST_SEGMENT_FRAMES + 1) * FIRST_SEGMENT_DT_NS
    tail = head[-1] + np.arange(1, n_frames - FIRST_SEGMENT_FRAMES + 1) * EXTENSION_DT_NS
    return np.concatenate([head, tail])


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
        raw = sorted(p for p in d.glob("*_analysis.dcd") if "pbcfix" not in p.name)
        fixed = sorted(d.glob("*_analysis_pbcfix.dcd"))
        top = sorted(d.glob("*analysis_topology.pdb"))
        if not (raw and top):
            continue
        import json
        js = sorted(d.glob("*_apo_rep*.json"))
        json_steps = 0
        if js:
            try:
                json_steps = int(json.loads(js[0].read_text()).get("md_production_steps_completed") or 0)
            except Exception:
                json_steps = 0
        steps = true_steps_completed(d, json_steps)
        if not steps:
            continue
        # Prefer the PBC-corrected sidecar; it is produced by the same path
        # the rest of the pipeline uses, so the imaging is identical.
        dcd = fixed[0] if fixed else raw[0]
        runs.append((d.parent.name, int(d.name.split("_")[-1]), dcd, top[0],
                     steps * 2e-6, bool(fixed)))   # 2 fs timestep -> ns
    if not runs:
        raise SystemExit("no usable apo runs found")

    rows = []
    for label, rep, dcd, top, total_ns, pre_imaged in runs:
        t = md.load(str(dcd), top=str(top))
        # Make the protein whole before measuring anything. Without this a
        # chain that drifts across the box boundary is recorded as scattered
        # fragments: in the G190E apo runs that inflated the pocket-lining
        # radius of gyration from 31 A to 79 A and the measured pocket volume
        # past 2000 A^3, which looks like a conformational change and is not.
        # A pbcfix sidecar has already had exactly this done to it.
        if not pre_imaged:
            try:
                t.make_molecules_whole(inplace=True)
                prot = t.topology.select("protein")
                if len(prot):
                    t.image_molecules(inplace=True,
                                      anchor_molecules=[{t.topology.atom(i) for i in prot}])
            except Exception as exc:
                print(f"  {label} rep{rep}: PBC correction failed ({exc})", flush=True)
        times = frame_times_ns(t.n_frames, total_ns)
        n_match = int((times <= args.match_ns).sum())
        if n_match == 0:
            n_match = t.n_frames
        u = mda.Universe(str(top))
        seg = _p66_segid(u)
        ca = u.select_atoms(f"segid {seg} and name CA and resid "
                            + " ".join(str(r + RESID_OFFSET) for r in NNIBP_P66))
        traj = t.xyz * 10.0                   # nm -> A

        for i in range(t.n_frames):
            u.atoms.positions = traj[i]
            window = "matched_10ns" if times[i] <= args.match_ns else "beyond_10ns"
            rows.append({
                "safe_label": label, "replicate": rep, "frame": i,
                "time_ns": float(times[i]), "window": window,
                "total_ns": total_ns,
                "pocket_volume": pocket_volume_proxy_from_universe(u, ""),
                "pocket_rg": float(np.sqrt((
                    (u.select_atoms(f"segid {seg} and name CA and resid "
                                    + " ".join(str(r + RESID_OFFSET) for r in NNIBP_P66)).positions
                     - ca.positions.mean(0)) ** 2).sum(-1).mean())),
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
        spacing = "uniform" if t.n_frames <= FIRST_SEGMENT_FRAMES else "two-stride"
        print(f"  {label:14s} rep{rep}  {t.n_frames:4d} frames  {total_ns:6.1f} ns  "
              f"{spacing:10s} matched={n_match} frames  RMSF={rmsf:.2f} A", flush=True)

    d = pd.DataFrame(rows)
    out = args.output_dir / "tables"
    out.mkdir(parents=True, exist_ok=True)

    # A replicate whose pocket has come apart cannot report a pocket volume.
    # G190E rep_02 reaches a pocket Ca Rg of 19.9 A and a median volume of
    # 2303 A^3, which on its own moved the G190E genotype mean from ~30 to
    # 868 A^3. Excluding these is a judgement about the simulation, not the
    # biology, so record it in the frame table and say so out loud.
    rep_rg = d.groupby(["safe_label", "replicate"]).pocket_rg.median()
    unfolded = set(rep_rg[rep_rg > POCKET_RG_MAX_A].index)
    d["pocket_intact"] = [
        (lbl, rep) not in unfolded for lbl, rep in zip(d.safe_label, d.replicate)
    ]
    d.to_csv(out / "apo_frame_metrics.csv", index=False)
    if unfolded:
        print(f"\n=== excluded: pocket Ca Rg > {POCKET_RG_MAX_A:.0f} A (pocket destroyed) ===")
        for lbl, rep in sorted(unfolded):
            print(f"  {lbl:14s} rep{rep}  pocket_rg={rep_rg[(lbl, rep)]:.1f} A  "
                  f"median volume={d[(d.safe_label == lbl) & (d.replicate == rep)].pocket_volume.median():.0f} A^3")

    d = d[d.pocket_intact]

    # A replicate shorter than the matched window has no matched window: every
    # frame it owns is inside it, so a 3.6 ns run would enter the table as a
    # 10 ns observation. Require the full window, and say which runs are too
    # short rather than quietly averaging them in.
    rep_total = d.groupby(["safe_label", "replicate"]).total_ns.first()
    too_short = sorted(rep_total[rep_total < args.match_ns * 0.98].index)
    if too_short:
        print(f"\n=== excluded from the matched window: shorter than {args.match_ns:.0f} ns ===")
        for lbl, rep in too_short:
            print(f"  {lbl:14s} rep{rep}  {rep_total[(lbl, rep)]:.1f} ns")
    long_enough = {k for k in rep_total.index if k not in set(too_short)}
    d["in_matched_window"] = [
        (lbl, rep) in long_enough for lbl, rep in zip(d.safe_label, d.replicate)
    ]

    m = d[(d.window == "matched_10ns") & d.in_matched_window]
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

    print(f"\n=== matched {args.match_ns:.0f} ns window, apo (intact replicates only) ===")
    print(summ.to_string())

    drift = []
    for lbl, g in d[d.total_ns > 20].groupby("safe_label"):
        head = g[g.window == "matched_10ns"].pocket_volume
        tail = g[g.window == "beyond_10ns"].pocket_volume
        if head.empty or tail.empty:
            continue
        a, b = head.mean(), tail.mean()
        drift.append({"safe_label": lbl, "n_rep": g.replicate.nunique(),
                      "first_10ns": round(a, 1),
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
