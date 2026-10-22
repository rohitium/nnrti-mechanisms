#!/usr/bin/env python3
"""Place the simulated apo ensembles between the crystallographic end states.

Reviewer 1 asks whether a cycle built on the DOR-bound crystal structure can
represent the unliganded enzyme at all, since the two differ and interconvert
slowly. That is testable rather than merely arguable: measure the NNIBP on
unliganded RT crystal structures with the same routine used for the
trajectories, and see where the apo simulations fall.

Residue identity is checked after renumbering, because the module's resids
carry a -3 offset relative to canonical HIV-1 RT numbering and a wrong offset
measures a different pocket instead of failing.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.compare_apo_to_crystal
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

#: canonical numbering -> expected residue, used to validate the offset
EXPECT = {97: "LEU", 100: "LYS", 103: "VAL", 178: "TYR", 185: "TYR", 187: "GLY"}
OFFSET = -3

STRUCTURES = [
    ("4NCG", "DOR-bound (holo)", "4NCG.pdb"),
    ("1DLO", "unliganded (apo)", "1DLO.pdb"),
    ("1HMV", "unliganded (apo)", "1HMV.pdb"),
]


def measure(path: Path, shift: int = OFFSET) -> tuple[float, int, list[str]]:
    import MDAnalysis as mda
    from nnrti.analysis.metrics import (pocket_volume_proxy_from_universe,
                                        _p66_segid, _NNBP_CA_RESIDS_P66)
    u = mda.Universe(str(path))
    if shift:
        u.residues.resids = u.residues.resids + shift
    seg = _p66_segid(u)
    bad = []
    for rid, name in EXPECT.items():
        sel = u.select_atoms(f"protein and name CA and segid {seg} and resid {rid}")
        got = sel.resnames[0] if sel.n_atoms else "MISSING"
        if got != name:
            bad.append(f"{rid}:{got}!={name}")
    n = u.select_atoms("protein and name CA and segid %s and (%s)" % (
        seg, " or ".join(f"resid {r}" for r in _NNBP_CA_RESIDS_P66))).n_atoms
    return float(pocket_volume_proxy_from_universe(u, "")), n, bad


def main(argv: list[str] | None = None) -> int:
    root = Path(__file__).resolve().parents[3]
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--structure-dir", type=Path,
                    default=root / "data/structures/apo_reference")
    ap.add_argument("--apo-summary", type=Path,
                    default=root / "results/analysis/apo_ensembles/tables/"
                                   "apo_matched_window_summary.csv")
    ap.add_argument("--output-dir", type=Path,
                    default=root / "results/analysis/apo_ensembles/tables")
    args = ap.parse_args(argv)

    rows = []
    for pdb_id, state, fname in STRUCTURES:
        p = args.structure_dir / fname
        if not p.exists():
            # *.pdb is gitignored, so these are fetched on demand rather than
            # committed. MDAnalysis cannot read the mmCIF form, so the PDB
            # format is the one that has to be present.
            import urllib.request
            url = f"https://files.rcsb.org/download/{fname}"
            print(f"fetching {pdb_id} from RCSB")
            p.parent.mkdir(parents=True, exist_ok=True)
            try:
                urllib.request.urlretrieve(url, p)
            except Exception as exc:
                print(f"  could not fetch {url}: {exc}")
                continue
        v, n, bad = measure(p)
        if bad:
            raise SystemExit(f"{pdb_id}: residue numbering mismatch {bad}")
        rows.append({"source": f"{pdb_id} crystal", "state": state,
                     "pocket_volume": round(v, 1), "n_ca": n})
        print(f"{pdb_id:6s} {state:20s} {n:2d}/15 CA   V = {v:6.1f} A^3")

    if args.apo_summary.exists():
        a = pd.read_csv(args.apo_summary)
        wt = a[a.safe_label == "wt"]
        if len(wt):
            r = wt.iloc[0]
            print(f"\nsimulated WT apo (matched 10 ns): "
                  f"{r.pocket_volume:.1f} +/- {r.pocket_volume_sd:.1f} A^3")
            rows.append({"source": "this work, WT apo ensemble",
                         "state": "unliganded (apo)",
                         "pocket_volume": round(float(r.pocket_volume), 1),
                         "n_ca": 15})

    args.output_dir.mkdir(parents=True, exist_ok=True)
    out = args.output_dir / "apo_crystal_reference.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"\nwrote {out.relative_to(root)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
