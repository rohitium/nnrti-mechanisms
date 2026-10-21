#!/usr/bin/env python3
"""DOR torsion distributions, WT-referenced.

Reviewer 2 (comment 3) asked for per-atom or torsion-level analyses to support
the claims about ligand displacement, which until now rested on contact counts
and centre-of-mass distances.

Four rotatable dihedrals define the DOR conformation:

    tau1  pyridinone-O ether bond
    tau2  O-chlorocyanobenzene bond
    tau3  pyridinone-N to the methylene linker
    tau4  methylene linker to the triazolinone

Atom names are those of the analysis topology, recovered from the torsion
definitions used before the September repository reorganisation.

These are intramolecular, so no periodic-image correction is needed: every atom
of the ligand stays in the same image as the rest of it regardless of where the
molecule sits in the box. That makes this the cleanest available probe of
whether a mutation changes the ligand's internal conformation as opposed to
moving it bodily.

Because torsions are circular, the summary uses circular means and circular
standard deviations; an arithmetic mean over angles that straddle +/-180 would
be meaningless.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.compute_dor_torsions
"""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

TORSIONS = [
    ("tau1", ["C12x", "C2x", "O1x", "C9x"]),
    ("tau2", ["C2x", "O1x", "C9x", "C10x"]),
    ("tau3", ["C4x", "N2x", "C15x", "C14x"]),
    ("tau4", ["N2x", "C15x", "C14x", "N5x"]),
]


def circ_mean_deg(a: np.ndarray) -> float:
    r = np.deg2rad(a[np.isfinite(a)])
    return float(np.rad2deg(np.arctan2(np.sin(r).mean(), np.cos(r).mean())))


def circ_std_deg(a: np.ndarray) -> float:
    r = np.deg2rad(a[np.isfinite(a)])
    R = np.hypot(np.sin(r).mean(), np.cos(r).mean())
    return float(np.rad2deg(np.sqrt(-2.0 * np.log(max(R, 1e-12)))))


def main(argv: list[str] | None = None) -> int:
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import calc_dihedrals

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--runs-root", type=Path, default=Path("results/md_runs"))
    p.add_argument("--output-dir", type=Path, default=Path("results/analysis/dor_torsions"))
    p.add_argument("--ligand", default="2KW")
    args = p.parse_args(argv)

    dirs = sorted(d for d in args.runs_root.glob("*/rep_*")
                  if "apo" not in d.parts and list(d.glob("*analysis.dcd")))
    rows = []
    for d in dirs:
        label, rep = d.parent.name, int(d.name.split("_")[-1])
        dcd = sorted(d.glob("*analysis.dcd"))[0]
        top = sorted(d.glob("*analysis_topology.pdb"))[0]
        try:
            u = mda.Universe(str(top), str(dcd))
        except Exception as e:
            print(f"  {label} rep{rep}: {e}", flush=True)
            continue
        lig = u.select_atoms(f"resname {args.ligand}")
        groups = {}
        ok = True
        for name, atoms in TORSIONS:
            sel = [lig.select_atoms(f"name {a}") for a in atoms]
            if any(s.n_atoms != 1 for s in sel):
                ok = False
                break
            groups[name] = sel
        if not ok:
            print(f"  {label} rep{rep}: torsion atoms not found, skipped", flush=True)
            continue
        for ts in u.trajectory:
            r = {"safe_label": label, "replicate": rep, "frame": ts.frame}
            for name, sel in groups.items():
                r[name] = float(np.degrees(calc_dihedrals(
                    sel[0].positions, sel[1].positions, sel[2].positions, sel[3].positions)[0]))
            rows.append(r)
        print(f"  {label:14s} rep{rep}  {len(u.trajectory)} frames", flush=True)

    if not rows:
        raise SystemExit("no torsions computed")
    d = pd.DataFrame(rows)
    out = args.output_dir / "tables"
    out.mkdir(parents=True, exist_ok=True)
    d.to_csv(out / "dor_torsions_per_frame.csv", index=False)

    recs = []
    for lbl, g in d.groupby("safe_label"):
        rec = {"safe_label": lbl, "n_frames": len(g)}
        for name, _ in TORSIONS:
            rec[f"{name}_mean"] = round(circ_mean_deg(g[name].to_numpy()), 1)
            rec[f"{name}_circsd"] = round(circ_std_deg(g[name].to_numpy()), 1)
        recs.append(rec)
    summ = pd.DataFrame(recs).set_index("safe_label").sort_index()
    summ.to_csv(out / "dor_torsions_summary.csv")

    print("\n=== circular mean +/- circular SD per torsion (degrees) ===")
    print(summ.to_string())
    if "wt" in summ.index:
        print("\n=== change in circular SD vs WT (degrees): broader = more flexible ===")
        dd = pd.DataFrame({n: summ[f"{n}_circsd"] - summ.loc["wt", f"{n}_circsd"]
                           for n, _ in TORSIONS}).drop(index="wt").round(1)
        print(dd.sort_values("tau4", ascending=False).to_string())
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
