#!/usr/bin/env python3
"""Place water in the NNIBP cavity left empty when DOR was deleted.

The apo systems were solvated independently of the holo ones, but the buried
NNIBP came out empty: one water lies within 6 A of the pocket centre, where a
~170 A^3 cavity should hold five or six. The protein then relaxes into a
vacuum, which is not a physical closing process and over-collapses the pocket
relative to the unliganded crystal structures. Ivetac & McCammon (JMB 2009,
388, 644) fill the vacated pocket with water for exactly this reason.

Rather than rebuild the topology, this relocates existing bulk waters into the
cavity. Each water is moved as a rigid body, so internal geometry and the
system's rigid-water constraints are untouched, the particle count does not
change, and the existing system XML still applies. Only coordinates change.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.fill_apo_cavity \
        --start results/md_runs/apo/wt/rep_01/assets/wt_apo_md_rep01_start.pdb \
        --out   filled_start.pdb
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

#: a water oxygen needs this much room from protein heavy atoms to sit without
#: clashing; roughly the O...C contact distance
MIN_PROTEIN_DIST = 3.2
#: and this much from another water oxygen
MIN_WATER_DIST = 2.8
#: bulk waters chosen for relocation must be at least this far from the protein,
#: so removing them leaves a hole that bulk fills in immediately
BULK_MIN_DIST = 12.0


def main(argv: list[str] | None = None) -> int:
    import openmm.app as app
    import openmm.unit as u
    from scipy.spatial import cKDTree
    from nnrti.analysis.metrics import (_NNBP_CA_RESIDS_P66, _NNBP_CA_RESID_P51)

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--start", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--radius", type=float, default=9.0,
                    help="search radius around the pocket centre, angstrom")
    ap.add_argument("--spacing", type=float, default=0.5)
    ap.add_argument("--max-waters", type=int, default=12)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args(argv)

    pdb = app.PDBFile(str(args.start))
    pos = np.asarray(pdb.positions.value_in_unit(u.angstrom))
    top = pdb.topology

    SOLV = {"HOH", "WAT", "TIP3", "SOL"}
    IONS = {"NA", "CL", "K", "MG", "SOD", "CLA"}
    prot_heavy, water_o, water_res = [], [], []
    for res in top.residues():
        nm = res.name.strip().upper()
        if nm in SOLV:
            idx = [a.index for a in res.atoms()]
            o = next((a.index for a in res.atoms()
                      if a.element is not None and a.element.symbol == "O"), None)
            if o is not None and len(idx) == 3:
                water_o.append(o)
                water_res.append(idx)
        elif nm not in IONS:
            prot_heavy += [a.index for a in res.atoms()
                           if a.element is not None and a.element.symbol != "H"]
    prot_heavy = np.asarray(prot_heavy, dtype=int)
    water_o = np.asarray(water_o, dtype=int)
    print(f"protein heavy atoms {prot_heavy.size}, waters {water_o.size}")

    # pocket centre: the Calpha cage used for every pocket volume in the paper
    # p66 and p51 share residue numbering, so the subunit has to be named
    # explicitly: matching numbers alone picks up both and puts the "pocket
    # centre" midway between the two subunits.
    want66 = set(_NNBP_CA_RESIDS_P66)
    prot_chain_size = {}
    for res in top.residues():
        if res.name.strip().upper() in SOLV | IONS:
            continue
        cid = res.chain.index
        prot_chain_size[cid] = prot_chain_size.get(cid, 0) + 1
    if len(prot_chain_size) < 2:
        raise SystemExit(f"expected two protein chains, found {len(prot_chain_size)}")
    ranked = sorted(prot_chain_size, key=prot_chain_size.get, reverse=True)
    p66, p51 = ranked[0], ranked[1]
    print(f"p66 = chain index {p66} ({prot_chain_size[p66]} residues), "
          f"p51 = chain index {p51} ({prot_chain_size[p51]})")
    ca = []
    for res in top.residues():
        if res.name.strip().upper() in SOLV | IONS:
            continue
        try:
            rid = int(res.id)
        except ValueError:
            continue
        hit = (res.chain.index == p66 and rid in want66) or \
              (res.chain.index == p51 and rid == _NNBP_CA_RESID_P51)
        if hit:
            ca += [a.index for a in res.atoms() if a.name == "CA"]
    ca = np.asarray(sorted(set(ca)), dtype=int)
    if not 14 <= ca.size <= 16:
        raise SystemExit(f"found {ca.size} pocket CA, expected 15 or 16; "
                         "check residue numbering and chain assignment")
    centre = pos[ca].mean(0)
    print(f"pocket centre from {ca.size} CA: {np.round(centre,1)}")

    # grid points with room for a water oxygen
    r, s = args.radius, args.spacing
    ax = [np.arange(centre[i] - r, centre[i] + r + s, s) for i in range(3)]
    grid = np.stack(np.meshgrid(*ax, indexing="ij"), -1).reshape(-1, 3)
    grid = grid[np.linalg.norm(grid - centre, axis=1) <= r]
    free = grid[cKDTree(pos[prot_heavy]).query(grid)[0] >= MIN_PROTEIN_DIST]
    if not len(free):
        raise SystemExit("no room for water in the cavity")
    free = free[cKDTree(pos[water_o]).query(free)[0] >= MIN_WATER_DIST]
    print(f"grid points with room: {len(free)}")

    # greedy, furthest-first so the waters spread through the cavity
    chosen = []
    order = np.argsort(np.linalg.norm(free - centre, axis=1))
    for i in order:
        p = free[i]
        if all(np.linalg.norm(p - q) >= MIN_WATER_DIST for q in chosen):
            chosen.append(p)
        if len(chosen) >= args.max_waters:
            break
    print(f"cavity sites selected: {len(chosen)}")
    if not chosen:
        print("nothing to do")
        return 0

    # donor waters from bulk, furthest from the protein
    d_prot = cKDTree(pos[prot_heavy]).query(pos[water_o])[0]
    cand = np.where(d_prot >= BULK_MIN_DIST)[0]
    rng = np.random.default_rng(args.seed)
    donors = rng.choice(cand, size=len(chosen), replace=False)

    new = pos.copy()
    for site, di in zip(chosen, donors):
        idx = water_res[di]
        shift = site - new[water_o[di]]
        new[idx] += shift                       # rigid translation
    print(f"relocated {len(chosen)} waters from bulk into the cavity")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        app.PDBFile.writeFile(top, new * u.angstrom, fh, keepIds=True)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
