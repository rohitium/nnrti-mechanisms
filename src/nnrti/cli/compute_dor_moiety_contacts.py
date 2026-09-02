#!/usr/bin/env python3
"""Protein heavy-atom contacts resolved by DOR moiety.

Why this exists
---------------
The manuscript reports a single number for the V106A dislocation -- "the total
number of protein heavy atoms within 4.5 A of DOR fell from 224 +/- 1 in WT to
(that figure predates the switch to a 4.0 A cutoff; see
manuscript/contact-cutoff-sensitivity.md)
212 +/- 1 in V106A". A whole-ligand count says the interface loosens but not
*where*, and DOR is not a homogeneous object: it is three ring systems joined by
an ether and a methylene, and the V106A slide toward Ser105 is directional, so
the loss should not be uniform across them.

This decomposes the same count, now at 4.0 A, over the ligand:

  chlorocyanophenyl   the ring that stacks against Tyr188
  pyridinone          the central ring bearing the Lys103 backbone H-bond
  triazolinone        the distal ring
  linker              ether O + methylene carbons not in any ring

and reports each as a mean +/- SEM over replicates, together with the
wild-type-referenced change. It also ranks RT residues by contact loss, so the
partner side of the disruption can be named rather than inferred.

Ring systems are identified from bond connectivity in the first frame, not from
atom names, because the MD topology renames the crystallographic ligand atoms.

Usage
-----
    PYTHONPATH=. python -m nnrti.cli.compute_dor_moiety_contacts
    PYTHONPATH=. python -m nnrti.cli.compute_dor_moiety_contacts \
        --mutations WT V106A --stride 5
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

LIG_RESNAME = "2KW"
CUTOFF_NM = 0.40  # 4.0 A (switched from 4.5 on 2026-08-31)
OFFSET = -3  # topology resid = canonical + OFFSET

V106A_SET = ["WT", "V106A", "V106A+F227L", "V106A+L234I", "V106A+P225H"]


def _remap(candidate: Path, repo_root: Path) -> Path:
    if candidate.exists():
        return candidate
    marker = "nnrti-mechanisms/"
    text = str(candidate)
    if marker in text:
        mapped = repo_root / text.split(marker, 1)[1]
        if mapped.exists():
            return mapped
    return candidate


def _moieties(traj) -> dict[str, list[int]]:
    """Partition DOR heavy atoms into ring systems + linker, from connectivity."""
    import networkx as nx

    top = traj.topology
    lig = [a for a in top.atoms if a.residue.name == LIG_RESNAME and a.element.symbol != "H"]
    if not lig:
        raise ValueError("no ligand heavy atoms found")
    idx = {a.name: a.index for a in lig}
    xyz = traj.xyz[0] * 10.0  # nm -> A

    G = nx.Graph()
    G.add_nodes_from(idx)
    for i, j in itertools.combinations(idx, 2):
        if np.linalg.norm(xyz[idx[i]] - xyz[idx[j]]) < 1.85:
            G.add_edge(i, j)
    rings = [sorted(c) for c in nx.cycle_basis(G) if len(c) >= 5]

    out: dict[str, list[int]] = {}
    for r in rings:
        nbrs = {n for a in r for n in G[a] if n not in r}
        # Every ring takes its own exocyclic substituents (Cl and nitrile on the
        # phenyl; CF3 and carbonyl on the pyridinone; methyl and carbonyl on the
        # triazolinone). Treating only the phenyl this way -- as an earlier
        # version did -- understated the pyridinone and triazolinone surfaces and
        # made the two scripts disagree. The shared ether oxygen (O1) is left to
        # the linker so it is not double counted.
        name = ("chlorocyanophenyl" if any(n.startswith("Cl") for n in nbrs)
                else ("triazolinone" if len(r) == 5 else "pyridinone"))
        sub = list(r) + [n for n in nbrs
                         if len(list(G[n])) <= 2 and not n.startswith("O1")]
        out[name] = sorted({idx[a] for a in sub})

    assigned = {i for v in out.values() for i in v}
    out["linker"] = sorted(set(idx.values()) - assigned)
    out["whole_ligand"] = sorted(idx.values())
    return out


def _protein_heavy(traj) -> np.ndarray:
    top = traj.topology
    return np.array(
        [
            a.index
            for a in top.atoms
            if a.residue.name != LIG_RESNAME
            and a.element.symbol not in ("H",)
            and a.residue.is_protein
        ],
        dtype=int,
    )


def _counts(traj, lig_idx: np.ndarray, prot_idx: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Per-frame contact counts under both conventions.

    Returns (pair_count, atom_count):

      pair_count  number of (protein heavy atom, ligand heavy atom) pairs closer
                  than the cutoff. This is what `_ncontacts` in
                  compute_mechanism_coordinates.py reports, and therefore what
                  the manuscript's existing "224" and "45.9" figures are, despite
                  both being described there as counts of atoms.
      atom_count  number of *distinct* protein heavy atoms within the cutoff of
                  any ligand atom -- the quantity the manuscript's wording
                  actually describes.

    Both are reported so the two can never be silently conflated again.
    """
    import mdtraj as md

    pairs = np.array(list(itertools.product(prot_idx, lig_idx)), dtype=int)
    d = md.compute_distances(traj, pairs, periodic=True)
    d = d.reshape(traj.n_frames, len(prot_idx), len(lig_idx))
    close = d < CUTOFF_NM
    return (close.sum(axis=(1, 2)).astype(float),
            close.any(axis=2).sum(axis=1).astype(float))


def _per_residue(traj, lig_idx: np.ndarray, prot_idx: np.ndarray) -> dict[str, float]:
    """Mean per-frame contacting-heavy-atom count, keyed by canonical residue."""
    import mdtraj as md

    top = traj.topology
    pairs = np.array(list(itertools.product(prot_idx, lig_idx)), dtype=int)
    d = md.compute_distances(traj, pairs, periodic=True)
    d = d.reshape(traj.n_frames, len(prot_idx), len(lig_idx))
    hit = (d.min(axis=2) < CUTOFF_NM)  # frames x prot_atoms
    per_atom = hit.mean(axis=0)
    out: dict[str, float] = {}
    for a_i, val in zip(prot_idx, per_atom):
        if val <= 0:
            continue
        res = top.atom(int(a_i)).residue
        key = f"{res.name}{res.resSeq - OFFSET}"
        out[key] = out.get(key, 0.0) + float(val)
    return out


def main() -> int:
    import mdtraj as md

    repo_root = Path(__file__).resolve().parents[3]
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--manifest", type=Path, default=repo_root / "manifests/md_manifest.csv")
    ap.add_argument("--mutations", nargs="*", default=V106A_SET)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--output-dir", type=Path,
                    default=repo_root / "results/analysis/mechanisms")
    args = ap.parse_args()

    man = pd.read_csv(args.manifest)
    man = man[man.mutation.isin(args.mutations)]

    rows, resrows = [], []
    for _, r in man.sort_values(["mutation", "replicate"]).iterrows():
        js = _remap(Path(str(r["output_json"])), repo_root)
        data = json.loads(js.read_text())
        topo = _remap(Path(str(data.get("analysis_topology_pdb") or "").strip()), repo_root)
        dcd = _remap(Path(str(data.get("analysis_dcd") or "").strip()), repo_root)
        if not (topo.exists() and dcd.exists()):
            print(f"  skip {r['mutation']} rep{r['replicate']}: missing trajectory")
            continue

        traj = md.load(str(dcd), top=str(topo), stride=args.stride)
        traj.make_molecules_whole(inplace=True)  # PBC: bond-graph traversal
        moi = _moieties(traj)
        prot = _protein_heavy(traj)

        rec = {"mutation": r["mutation"], "replicate": int(r["replicate"]),
               "n_frames": traj.n_frames}
        for name, atoms in moi.items():
            pc, ac = _counts(traj, np.asarray(atoms, dtype=int), prot)
            rec[name] = float(pc.mean())
            rec[f"{name}_sd"] = float(pc.std(ddof=1))
            rec[f"{name}_atoms"] = float(ac.mean())
        rows.append(rec)
        print(f"  {r['mutation']:14s} rep{r['replicate']}  "
              f"whole {rec['whole_ligand']:6.1f}  "
              f"chl {rec['chlorocyanophenyl']:5.1f}  "
              f"pyr {rec['pyridinone']:5.1f}  "
              f"tri {rec['triazolinone']:5.1f}  "
              f"lnk {rec['linker']:5.1f}")

        for res, val in _per_residue(traj, np.asarray(moi["whole_ligand"], int), prot).items():
            resrows.append({"mutation": r["mutation"], "replicate": int(r["replicate"]),
                            "residue": res, "contacts": val})

    df = pd.DataFrame(rows)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_dir / "dor_moiety_contacts_per_replicate.csv", index=False)

    moieties = ["whole_ligand", "chlorocyanophenyl", "pyridinone", "triazolinone", "linker"]
    cols = moieties + [f"{c}_atoms" for c in moieties]
    agg = df.groupby("mutation")[cols].agg(["mean", "sem"])
    wt = agg.loc["WT"] if "WT" in agg.index else None

    summary = []
    for m in agg.index:
        rec = {"mutation": m}
        for c in cols:
            rec[c] = round(float(agg.loc[m, (c, "mean")]), 1)
            rec[f"{c}_sem"] = round(float(agg.loc[m, (c, "sem")]), 1)
            if wt is not None:
                rec[f"d_{c}"] = round(float(agg.loc[m, (c, "mean")] - wt[(c, "mean")]), 1)
        summary.append(rec)
    sdf = pd.DataFrame(summary)
    sdf.to_csv(args.output_dir / "dor_moiety_contacts_summary.csv", index=False)

    print("\n=== ATOM-PAIR contacts within 4.0 A (the manuscript's convention) ===")
    hdr = f"{'genotype':16s}" + "".join(f"{c[:12]:>18s}" for c in moieties)
    print(hdr)
    for _, rec in sdf.iterrows():
        line = f"{rec['mutation']:16s}"
        for c in moieties:
            line += f"{rec[c]:8.1f} ± {rec[f'{c}_sem']:<4.1f}   "[:18]
        print(line)
    print("\n=== DISTINCT protein heavy atoms within 4.0 A ===")
    print(f"{'genotype':16s}" + "".join(f"{c[:12]:>18s}" for c in moieties))
    for _, rec in sdf.iterrows():
        line = f"{rec['mutation']:16s}"
        for c in moieties:
            line += f"{rec[c + '_atoms']:8.1f} ± {rec[f'{c}_atoms_sem']:<4.1f}   "[:18]
        print(line)

    if wt is not None:
        for label, suffix in (("ATOM-PAIR contacts", ""), ("DISTINCT atoms", "_atoms")):
            print(f"\n=== change vs WT, {label} ===")
            print(f"{'genotype':16s}" + "".join(f"{c[:12]:>14s}" for c in moieties))
            for _, rec in sdf.iterrows():
                if rec["mutation"] == "WT":
                    continue
                print(f"{rec['mutation']:16s}"
                      + "".join(f"{rec[f'd_{c}{suffix}']:+13.1f} " for c in moieties))

    if resrows:
        rdf = pd.DataFrame(resrows)
        rep = rdf.groupby(["mutation", "residue"]).contacts.mean().unstack(0).fillna(0.0)
        rdf.to_csv(args.output_dir / "dor_residue_contacts_per_replicate.csv", index=False)
        if "WT" in rep.columns:
            v = [c for c in rep.columns if c != "WT"]
            delta = rep[v].mean(axis=1) - rep["WT"]
            print("\n=== RT residues losing the most contact with DOR "
                  "(mean over V106A genotypes vs WT) ===")
            for res, d in delta.sort_values().head(10).items():
                print(f"  {res:10s} {rep.loc[res,'WT']:6.2f} -> "
                      f"{rep.loc[res, v].mean():6.2f}   ({d:+.2f})")
            print("\n=== residues gaining the most ===")
            for res, d in delta.sort_values(ascending=False).head(5).items():
                print(f"  {res:10s} {rep.loc[res,'WT']:6.2f} -> "
                      f"{rep.loc[res, v].mean():6.2f}   ({d:+.2f})")
            delta.sort_values().to_csv(args.output_dir / "dor_residue_contact_delta.csv")

    print(f"\nWrote {args.output_dir}/dor_moiety_contacts_summary.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
