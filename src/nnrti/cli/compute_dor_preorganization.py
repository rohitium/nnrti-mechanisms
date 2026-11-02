#!/usr/bin/env python3
"""What DOR pays, on its own, to adopt its bound conformation.

Reviewer 1 asks whether mutation-induced remodelling is poorly reflected in
ddG_bind. One term ddG_bind omits by construction is the cost of pre-organising
the ligand: single-trajectory MM/GBSA extracts the free ligand from the complex
frame, so the ligand is already in its bound conformation and pays nothing.

This measures that term the only way it can be measured reliably here -- by
counting, not by differencing energies. Absolute MM/GBSA energies are unusable
for this because a capped snapshot minimisation leaves the value dependent on
how far the minimiser descended (~2200 kcal/mol between 100 and 1000 iterations,
against a 7 kcal/mol panel range). A population ratio has no minimiser in it at
all:

    dG_preorg = -RT ln p,   p = P(free ligand is in the bound conformation)

The bound conformation is defined from the WT holo ensemble as a window of
+/- k circular SD about the circular mean of each of the four rotatable
dihedrals. Two corrections keep that definition honest:

  * capture: the same window is applied to the bound ensemble itself, and p is
    divided by the fraction it retains, so a window that clips the bound state
    does not inflate the cost.
  * conditioning: counting joint occupancy directly rests on a handful of
    frames. The four torsions are close to independent (verified, not assumed),
    so the product of the four per-torsion fractions -- each resting on
    thousands of frames -- is reported alongside it as the better-conditioned
    estimate. The two are required to agree.

The quantity is a property of the free ligand, so it is one number for the whole
panel, not a per-genotype term: free DOR does not know which mutant it will
bind. It therefore sets the size of a term ddG_bind omits; it does not explain
differences between genotypes.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.compute_dor_preorganization
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

#: kcal/mol at the production temperature.
RT_KCAL = 0.001987204 * 300.0

TAUS = ["tau1", "tau2", "tau3", "tau4"]

#: Window half-widths, in units of the bound-state circular SD. Reported as a
#: sweep rather than a single choice because the free-ligand density keeps
#: rising as the window widens -- the edge of the bound basin is genuinely
#: fuzzy, and that spread is the dominant uncertainty, not counting noise.
K_VALUES = (1.5, 2.0, 2.5, 3.0)

#: Frame spacing of the solvent runs: 2000 analysis frames over 100 ns.
FRAME_DT_NS = 0.05

#: A rotamer change must persist this many frames to count as a crossing, so
#: thermal rattling at a barrier top is not mistaken for a transition.
MIN_WELL_FRAMES = 5

#: Below this many crossings in a trajectory, a torsion's population carries the
#: memory of the starting conformer and the joint estimate inherits that bias.
#: This is not a tunable: it is the point at which the number of independent
#: samples drops into single figures.
MIN_CROSSINGS = 10


def circ_mean_deg(a: np.ndarray) -> float:
    r = np.deg2rad(np.asarray(a, float))
    r = r[np.isfinite(r)]
    return float(np.rad2deg(np.arctan2(np.sin(r).mean(), np.cos(r).mean())))


def circ_sd_deg(a: np.ndarray) -> float:
    r = np.deg2rad(np.asarray(a, float))
    r = r[np.isfinite(r)]
    R = np.hypot(np.sin(r).mean(), np.cos(r).mean())
    return float(np.rad2deg(np.sqrt(-2.0 * np.log(max(R, 1e-12)))))


def circ_dist_deg(a: np.ndarray, centre: float) -> np.ndarray:
    """Shortest angular separation, so a window never straddles the +/-180 seam."""
    return np.abs((np.asarray(a, float) - centre + 180.0) % 360.0 - 180.0)



def count_crossings(angles_deg: np.ndarray) -> tuple[int, np.ndarray]:
    """Rotamer crossings and well occupancies for one torsion of one replicate.

    Wells are 120 degree sectors. Runs shorter than MIN_WELL_FRAMES are dropped
    before counting, so a torsion sitting on a barrier top does not register a
    crossing on every frame.
    """
    a = ((np.asarray(angles_deg, float) + 180.0) % 360.0) - 180.0
    well = (np.floor((a + 180.0) / 120.0).astype(int)) % 3
    edges = np.concatenate(([0], np.flatnonzero(np.diff(well)) + 1, [len(well)]))
    segs = [(well[edges[i]], edges[i + 1] - edges[i]) for i in range(len(edges) - 1)]
    kept = [w for w, n in segs if n >= MIN_WELL_FRAMES]
    crossings = sum(1 for i in range(1, len(kept)) if kept[i] != kept[i - 1])
    occupancy = np.bincount(well, minlength=3) / len(well)
    return crossings, occupancy


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--torsions-csv", type=Path,
                   default=Path("results/analysis/dor_torsions_solution/tables/"
                                "dor_torsions_per_frame.csv"))
    p.add_argument("--solution-label", default="dor_solvent")
    p.add_argument("--bound-label", default="wt")
    p.add_argument("--output-dir", type=Path,
                   default=Path("results/analysis/dor_preorganization"))
    args = p.parse_args(argv)

    d = pd.read_csv(args.torsions_csv)
    sol = d[d.safe_label == args.solution_label]
    bound = d[d.safe_label == args.bound_label]
    if sol.empty or bound.empty:
        raise SystemExit(
            f"need both {args.solution_label!r} and {args.bound_label!r} in "
            f"{args.torsions_csv}; found {sorted(d.safe_label.unique())}")

    ref = {t: (circ_mean_deg(bound[t]), circ_sd_deg(bound[t])) for t in TAUS}

    print(f"free ligand : {len(sol):5d} frames "
          f"({sol.replicate.nunique()} replicates)")
    print(f"WT bound    : {len(bound):5d} frames "
          f"({bound.replicate.nunique()} replicates)\n")

    print("=== binding rigidifies every torsion ===")
    print(f"{'':6s}{'bound mean':>12s}{'bound SD':>10s}{'free SD':>9s}{'ratio':>8s}")
    for t in TAUS:
        m, s = ref[t]
        fs = circ_sd_deg(sol[t])
        print(f"{t:6s}{m:12.1f}{s:10.1f}{fs:9.1f}{fs / s:7.1f}x")

    rows = []
    for k in K_VALUES:
        per_tau, joint_sol = {}, np.ones(len(sol), bool)
        for t in TAUS:
            m, s = ref[t]
            ok = circ_dist_deg(sol[t].to_numpy(), m) <= k * s
            per_tau[t] = float(ok.mean())
            joint_sol &= ok

        joint_bound = np.ones(len(bound), bool)
        for t in TAUS:
            m, s = ref[t]
            joint_bound &= circ_dist_deg(bound[t].to_numpy(), m) <= k * s
        capture = float(joint_bound.mean())

        n_hit = int(joint_sol.sum())
        p_direct = n_hit / len(sol)
        p_indep = float(np.prod([per_tau[t] for t in TAUS]))
        # Capture-correct the better-conditioned estimate.
        p_corr = p_indep / capture if capture > 0 else np.nan
        dg = -RT_KCAL * np.log(p_corr) if p_corr > 0 else np.nan

        rows.append({
            "k": k, "n_hits": n_hit, "n_solution_frames": len(sol),
            "p_direct": p_direct, "p_independent": p_indep,
            "bound_capture": capture, "p_corrected": p_corr,
            "dG_preorg_kcal": dg,
            **{f"p_{t}": per_tau[t] for t in TAUS},
        })

    res = pd.DataFrame(rows)

    print("\n=== fraction of the free ensemble in the bound conformation ===")
    print(f"{'k':>5s}" + "".join(f"{t:>8s}" for t in TAUS)
          + f"{'direct':>9s}{'indep':>9s}{'capture':>9s}{'corr':>9s}{'dG':>8s}")
    for r in rows:
        print(f"{r['k']:5.1f}"
              + "".join(f"{r[f'p_{t}'] * 100:7.1f}%" for t in TAUS)
              + f"{r['p_direct'] * 100:8.2f}%{r['p_independent'] * 100:8.2f}%"
              + f"{r['bound_capture'] * 100:8.1f}%{r['p_corrected'] * 100:8.2f}%"
              + f"{r['dG_preorg_kcal']:8.2f}")

    # The product estimate is only legitimate if the torsions really are close
    # to independent. State the discrepancy rather than trusting it silently.
    worst = max(abs(r["p_independent"] / r["p_direct"] - 1.0)
                for r in rows if r["p_direct"] > 0)
    print(f"\nindependence check: product vs direct count differs by at most "
          f"{worst * 100:.0f}% across the sweep")
    if worst > 1.0:
        print("  WARNING: torsions are not independent; trust the direct count")

    lo, hi = res.dG_preorg_kcal.min(), res.dG_preorg_kcal.max()
    print(f"\ndG_preorg = {lo:.1f} to {hi:.1f} kcal/mol "
          f"(window-width spread, the dominant uncertainty)")

    print("\n=== per-replicate reproducibility (k=2.0) ===")
    for rep, g in sol.groupby("replicate"):
        j = np.ones(len(g), bool)
        for t in TAUS:
            m, s = ref[t]
            j &= circ_dist_deg(g[t].to_numpy(), m) <= 2.0 * s
        sds = " ".join(f"{circ_sd_deg(g[t]):3.0f}" for t in TAUS)
        print(f"  rep_{rep:02d}: {len(g):5d} frames  joint {j.mean() * 100:5.2f}%"
              f"  free SDs [{sds}]")

    # Convergence of the free ensemble. The joint population multiplies through
    # every torsion, so the worst-sampled one sets the real uncertainty on
    # dG_preorg -- which is why this runs every time rather than on request.
    print("\n=== convergence of the free ensemble: rotamer crossings in "
          f"{len(sol) * FRAME_DT_NS / sol.replicate.nunique():.0f} ns per replicate ===")
    diag = []
    for rep, g in sol.groupby("replicate"):
        g = g.sort_values("frame")
        for t in TAUS:
            n_cross, occ = count_crossings(g[t].to_numpy())
            ns_per = len(g) * FRAME_DT_NS / max(n_cross, 1)
            diag.append({"replicate": rep, "torsion": t, "crossings": n_cross,
                         "ns_per_crossing": ns_per,
                         **{f"well{i}_occupancy": occ[i] for i in range(3)}})
            flag = "  UNDERSAMPLED" if n_cross < MIN_CROSSINGS else ""
            print(f"  rep_{rep:02d} {t}: {n_cross:3d} crossings"
                  f"  ~{ns_per:5.1f} ns each"
                  f"  occupancy {occ[0]*100:5.1f}/{occ[1]*100:5.1f}/{occ[2]*100:5.1f}%"
                  f"{flag}")
    diag = pd.DataFrame(diag)

    worst = diag.loc[diag.crossings.idxmin()]
    if worst.crossings < MIN_CROSSINGS:
        print(f"\n  {worst.torsion} crosses only {int(worst.crossings)} times in the "
              f"worst replicate (~{worst.ns_per_crossing:.0f} ns per crossing).")
        print("  The joint population multiplies through it, so dG_preorg above is "
              "an order-of-magnitude")
        print("  estimate only -- quote it as 'a few kcal/mol', not to one decimal. "
              "Pinning it down")
        print("  needs enhanced sampling of the torsions, not more plain MD: at this "
              "crossing rate a")
        print("  converged population would take microseconds.")

    out = args.output_dir / "tables"
    out.mkdir(parents=True, exist_ok=True)
    res.to_csv(out / "dor_preorganization.csv", index=False)
    diag.to_csv(out / "dor_torsion_convergence.csv", index=False)
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
