#!/usr/bin/env python3
"""Sensitivity of the DOR contact analysis to the distance cutoff.

Why this exists
---------------
"Contact" has no single definition in the literature: 4.0 A and 4.5 A are both
widely used, and binding-site definitions often go to 5.0 A. Rather than defend
one choice, this recomputes the manuscript's contact quantities across a range
of cutoffs in a single pass over each trajectory, so the reader can see which
conclusions depend on the threshold and which do not.

Distances are computed once per replicate and thresholded at every cutoff, so
adding cutoffs is nearly free.

Usage
-----
    PYTHONPATH=. python -m src.analysis.cli.sweep_contact_cutoff
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

LIG = "2KW"
CUTOFFS_NM = (0.35, 0.40, 0.45, 0.50)
GENOTYPES = ["WT", "V106A", "V106A+F227L", "V106A+L234I", "V106A+P225H", "Y188L"]
KEY_RESIDUES = {103: "Lys103", 106: "Val106", 181: "Tyr181", 188: "Tyr188",
                190: "Gly190", 227: "Phe227", 229: "Trp229", 318: "Tyr318"}
OFFSET = -3


def _remap(p: Path, root: Path) -> Path:
    p = Path(str(p))
    if p.exists():
        return p
    if "nnrti-mechanisms/" in str(p):
        q = root / str(p).split("nnrti-mechanisms/", 1)[1]
        if q.exists():
            return q
    return p


def _moieties(traj) -> dict[str, list[int]]:
    import networkx as nx

    top = traj.topology
    lig = {a.name: a.index for a in top.atoms
           if a.residue.name == LIG and a.element.symbol != "H"}
    xyz = traj.xyz[0] * 10.0
    G = nx.Graph()
    G.add_nodes_from(lig)
    for i, j in itertools.combinations(lig, 2):
        if np.linalg.norm(xyz[lig[i]] - xyz[lig[j]]) < 1.85:
            G.add_edge(i, j)
    out: dict[str, list[int]] = {}
    for ring in [sorted(c) for c in nx.cycle_basis(G) if len(c) >= 5]:
        nb = {n for a in ring for n in G[a] if n not in ring}
        nm = ("chlorocyanophenyl" if any(n.startswith("Cl") for n in nb)
              else ("triazolinone" if len(ring) == 5 else "pyridinone"))
        ext = list(ring) + [n for n in nb
                            if len(list(G[n])) <= 2 and not n.startswith("O1")]
        out[nm] = sorted({lig[a] for a in ext})
    assigned = {i for v in out.values() for i in v}
    out["linker"] = sorted(set(lig.values()) - assigned)
    out["whole_ligand"] = sorted(lig.values())
    return out


def main() -> int:
    import mdtraj as md

    root = Path(__file__).resolve().parents[3]
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--manifest", type=Path, default=root / "manifests/md_manifest.csv")
    ap.add_argument("--mutations", nargs="*", default=GENOTYPES)
    ap.add_argument("--stride", type=int, default=10)
    ap.add_argument("--output-dir", type=Path, default=root / "results/analysis/mechanisms")
    args = ap.parse_args()

    man = pd.read_csv(args.manifest)
    man = man[man.mutation.isin(args.mutations)]

    rows, resrows = [], []
    for _, r in man.sort_values(["mutation", "replicate"]).iterrows():
        js = _remap(Path(str(r["output_json"])), root)
        data = json.loads(js.read_text())
        topo = _remap(Path(str(data.get("analysis_topology_pdb") or "")), root)
        dcd = _remap(Path(str(data.get("analysis_dcd") or "")), root)
        if not (topo.exists() and dcd.exists()):
            continue
        traj = md.load(str(dcd), top=str(topo), stride=args.stride)
        traj.make_molecules_whole(inplace=True)
        top = traj.topology
        moi = _moieties(traj)
        prot = np.array([a.index for a in top.atoms
                         if a.residue.name != LIG and a.element.symbol != "H"], dtype=int)

        for name, atoms in moi.items():
            la = np.asarray(atoms, dtype=int)
            pairs = np.array(list(itertools.product(prot, la)), dtype=int)
            d = md.compute_distances(traj, pairs, periodic=True)
            d = d.reshape(traj.n_frames, len(prot), len(la))
            for cut in CUTOFFS_NM:
                close = d < cut
                rows.append(dict(mutation=r["mutation"], replicate=int(r["replicate"]),
                                 moiety=name, cutoff_A=round(cut * 10, 1),
                                 pairs=float(close.sum(axis=(1, 2)).mean()),
                                 atoms=float(close.any(axis=2).sum(axis=1).mean())))

        if r["mutation"] in ("WT",):
            for can, nm in KEY_RESIDUES.items():
                pa = np.array([a.index for a in top.atoms
                               if a.residue.resSeq == can + OFFSET
                               and a.residue.name != LIG
                               and a.element.symbol != "H"
                               and a.residue.chain.index == 0], dtype=int)
                if pa.size == 0:
                    continue
                la = np.asarray(moi["whole_ligand"], dtype=int)
                pairs = np.array(list(itertools.product(pa, la)), dtype=int)
                d = md.compute_distances(traj, pairs, periodic=True)
                for cut in CUTOFFS_NM:
                    resrows.append(dict(residue=nm, replicate=int(r["replicate"]),
                                        cutoff_A=round(cut * 10, 1),
                                        pairs=float((d < cut).sum(axis=1).mean())))
        print(f"  {r['mutation']:14s} rep{r['replicate']} done")

    df = pd.DataFrame(rows)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_dir / "contact_cutoff_sweep.csv", index=False)

    agg = df.groupby(["mutation", "moiety", "cutoff_A"]).pairs.mean().reset_index()
    piv = agg.pivot_table(index=["moiety", "mutation"], columns="cutoff_A", values="pairs")

    print("\n=== atom-pair contacts vs cutoff ===")
    for moi in ["whole_ligand", "chlorocyanophenyl", "pyridinone", "triazolinone"]:
        if moi not in piv.index.get_level_values(0):
            continue
        print(f"\n-- {moi} --")
        sub = piv.loc[moi]
        print(sub.round(1).to_string())
        if "WT" in sub.index:
            print("   change vs WT (%):")
            for m in sub.index:
                if m == "WT":
                    continue
                pct = 100 * (sub.loc[m] - sub.loc["WT"]) / sub.loc["WT"]
                print(f"     {m:14s}" + "".join(f"{v:+8.1f}%" for v in pct))

    if resrows:
        rd = pd.DataFrame(resrows).groupby(["residue", "cutoff_A"]).pairs.mean().unstack()
        rd.to_csv(args.output_dir / "contact_cutoff_sweep_residues.csv")
        print("\n=== WT: key residue contacts with whole ligand, vs cutoff ===")
        print(rd.round(1).to_string())

    print(f"\nWrote {args.output_dir}/contact_cutoff_sweep.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
