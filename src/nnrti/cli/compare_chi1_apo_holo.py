#!/usr/bin/env python3
"""Does a mutation move Tyr181/Tyr188 in the unliganded enzyme as well as the bound one?

Reviewer 1 asks (comment 3) whether mutation-induced changes occur in both apo
and holo RT, and whether such remodeling is therefore poorly reflected in
ddG_bind. The thermodynamic cycle introduces each mutation twice, once in each
state, so a mutation that perturbs both states in the same direction cancels
and ddG_bind cannot see it. That is a measurement, not a hypothesis, once the
same observable is available in both states:

    d_holo = X(mutant, holo) - X(wild type, holo)
    d_apo  = X(mutant, apo)  - X(wild type, apo)

A genotype on the diagonal d_apo = d_holo perturbs both states equally and is
invisible to the cycle. One off the diagonal has a differential effect the
cycle should register.

The observable is the chi1 dihedral of Tyr181 and Tyr188. The reviewer names
that rotation as what opening the cleft involves, and unlike a pocket volume it
is a property of the side chain alone -- no cavity definition, no solvation
assumption, and no sensitivity to how the apo cavity was equilibrated.

chi1 is circular, so means, spreads and differences all use circular
statistics; an ordinary mean of angles near the +/-180 wrap is meaningless.
A torsion is intra-residue, so no periodic-image correction is needed.

Genotypes that remove the residue are excluded from its panel: Y181C has no
Tyr181, Y188L no Tyr188.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.compare_chi1_apo_holo
"""
from __future__ import annotations

import argparse
import glob
import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[3]

RESIDUES = {"Tyr181": 181, "Tyr188": 188}
RESID_OFFSET = -3
CHI1 = ("N", "CA", "CB", "CG")
MISMATCH = {"k103n_m230l", "a98g_f227c", "g190e"}


def circular_mean(deg: np.ndarray) -> float:
    r = np.deg2rad(deg)
    return float(np.rad2deg(np.arctan2(np.sin(r).mean(), np.cos(r).mean())))


def circular_sd(deg: np.ndarray) -> float:
    """Circular standard deviation in degrees (Mardia)."""
    r = np.deg2rad(deg)
    R = np.hypot(np.sin(r).mean(), np.cos(r).mean())
    R = min(max(R, 1e-12), 1.0)
    return float(np.rad2deg(np.sqrt(-2.0 * np.log(R))))


def circular_delta(a: float, b: float) -> float:
    """Signed smallest difference a - b, on (-180, 180]."""
    return float((a - b + 180.0) % 360.0 - 180.0)


def dihedral(p: np.ndarray) -> float:
    b0, b1, b2 = p[0] - p[1], p[2] - p[1], p[3] - p[2]
    b1 /= np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))


def chi1_series(dcd: str, topo: str, resid: int) -> np.ndarray | None:
    """chi1 per frame for one residue of p66, or None if the side chain is absent."""
    import MDAnalysis as mda

    from nnrti.analysis.metrics import _p66_segid

    u = mda.Universe(topo, dcd)
    seg = _p66_segid(u)
    sel = u.select_atoms(f"segid {seg} and resid {resid + RESID_OFFSET}")
    if sel.n_atoms == 0:
        return None
    # The residue must still be a tyrosine. Leucine also has N-CA-CB-CG, so an
    # atom-presence test alone lets Y188L through and silently compares a
    # leucine chi1 against a tyrosine's.
    resnames = set(sel.resnames)
    if resnames != {"TYR"}:
        return None
    names = set(sel.names)
    if not {"N", "CA", "CB", "CG"} <= names:
        return None
    groups = [sel.select_atoms(f"name {n}") for n in CHI1]
    if any(g.n_atoms != 1 for g in groups):
        return None
    out = np.empty(u.trajectory.n_frames)
    for i, _ in enumerate(u.trajectory):
        out[i] = dihedral(np.array([g.positions[0] for g in groups]))
    return out


def replicate_dirs(state: str, genotype: str) -> list[str]:
    if state == "apo":
        return sorted(glob.glob(f"results/md_runs/apo/{genotype}/rep_*"))
    # holo directories keep the manuscript's capitalisation
    for name in (genotype.upper().replace("_", "+"), genotype, "wt"):
        hits = sorted(glob.glob(f"results/md_runs/{name}/rep_*"))
        if hits:
            return hits
    return []


def trajectory_in(d: str) -> tuple[str, str] | None:
    fix = glob.glob(os.path.join(d, "*_analysis_pbcfix.dcd"))
    raw = [p for p in glob.glob(os.path.join(d, "*_analysis.dcd")) if "pbcfix" not in p]
    topo = [p for p in glob.glob(os.path.join(d, "*_analysis_topology.pdb"))
            if "from_md" not in os.path.basename(p)]
    src = (fix or raw)
    if not src or not topo:
        return None
    return sorted(src)[-1], sorted(topo)[-1]


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--output-dir", type=Path,
                    default=REPO_ROOT / "results/analysis/chi1_apo_holo")
    ap.add_argument("--genotypes", nargs="*", default=None)
    args = ap.parse_args(argv)
    tables = args.output_dir / "tables"
    tables.mkdir(parents=True, exist_ok=True)

    genotypes = args.genotypes or sorted(
        os.path.basename(os.path.dirname(d))
        for d in glob.glob("results/md_runs/apo/*/rep_01")
    )

    rows = []
    for genotype in genotypes:
        for state in ("apo", "holo"):
            for d in replicate_dirs(state, genotype):
                found = trajectory_in(d)
                if not found:
                    continue
                dcd, topo = found
                rep = int(os.path.basename(d).split("_")[-1])
                for label, resid in RESIDUES.items():
                    series = chi1_series(dcd, topo, resid)
                    if series is None:
                        continue
                    rows.append(dict(
                        genotype=genotype, state=state, replicate=rep, residue=label,
                        n_frames=len(series),
                        chi1_mean=round(circular_mean(series), 2),
                        chi1_sd=round(circular_sd(series), 2),
                    ))
                logging.info("%s %s rep%02d done", genotype, state, rep)

    per_rep = pd.DataFrame(rows)
    per_rep.to_csv(tables / "chi1_per_replicate.csv", index=False)

    # genotype-level: circular mean of the replicate means, unit = replicate
    agg = []
    for (g, s, r), grp in per_rep.groupby(["genotype", "state", "residue"]):
        agg.append(dict(genotype=g, state=s, residue=r, n_rep=len(grp),
                        chi1=round(circular_mean(grp.chi1_mean.to_numpy()), 2),
                        spread=round(circular_sd(grp.chi1_mean.to_numpy()), 2)))
    genotype_level = pd.DataFrame(agg)
    genotype_level.to_csv(tables / "chi1_per_genotype.csv", index=False)

    # the comparison the reviewer asks for
    wide = genotype_level.pivot_table(index=["genotype", "residue"], columns="state",
                                      values="chi1", aggfunc="first")
    out = []
    for residue in RESIDUES:
        try:
            wt_apo = wide.loc[("wt", residue), "apo"]
            wt_holo = wide.loc[("wt", residue), "holo"]
        except KeyError:
            continue
        for (g, r), row in wide.iterrows():
            if r != residue or g == "wt":
                continue
            if pd.isna(row.get("apo")) or pd.isna(row.get("holo")):
                continue
            out.append(dict(
                genotype=g, residue=residue,
                delta_apo=round(circular_delta(row["apo"], wt_apo), 2),
                delta_holo=round(circular_delta(row["holo"], wt_holo), 2),
                mismatch_genotype=g in MISMATCH,
            ))
    comparison = pd.DataFrame(out)
    comparison["delta_difference"] = (comparison.delta_apo - comparison.delta_holo).round(2)
    comparison.to_csv(tables / "chi1_apo_vs_holo.csv", index=False)

    for residue in RESIDUES:
        sub = comparison[comparison.residue == residue]
        if sub.empty:
            continue
        print(f"\n=== {residue}: shift from wild type, degrees ===")
        print(sub.sort_values("delta_holo")[
            ["genotype", "delta_apo", "delta_holo", "delta_difference", "mismatch_genotype"]
        ].to_string(index=False))
    logging.info("wrote %s", tables)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
