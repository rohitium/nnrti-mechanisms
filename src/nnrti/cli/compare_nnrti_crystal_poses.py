#!/usr/bin/env python3
"""Where Tyr181 and Tyr188 sit in RT bound to different NNRTIs.

The manuscript argues that doravirine's aromatic packing rests on Tyr188 alone,
because the Tyr181 side chain is rotated away from the drug. That claim is made
against the crystallographic observation for 4NCG. This compares the two
tyrosines across four NNRTI-bound structures of wild-type RT, so the
Introduction can say how DOR differs from the first- and second-generation
drugs rather than asserting it:

    4NCG   2.58 A X-ray     doravirine    (2KW)
    7Z2G   3.65 A cryo-EM   doravirine    (2KW) + 38-mer DNA
    2ZD1   1.80 A X-ray     rilpivirine   (T27)
    1FK9   2.50 A X-ray     efavirenz     (EFZ)
    1VRT   2.20 A X-ray     nevirapine    (NVP)

Everything is reported after superposition on the 16 NNIBP-lining Calpha atoms
of the reference, so the numbers describe the side chains inside a common
pocket frame and not a whole-enzyme difference. The ligand-contact distances
are computed with no superposition at all.

Run:
    PYTHONPATH=src python -m nnrti.cli.compare_nnrti_crystal_poses
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from nnrti.cli.compare_dor_crystal_poses import (
    NNIBP_P51,
    NNIBP_P66,
    Struct,
    apply,
    dihedral,
    ensure_pdb,
    kabsch,
    pocket_volume,
    rmsd,
)

# label -> (pdb id, ligand component, method, resolution A, nucleic acid)
STRUCTURES: dict[str, tuple[str, str, str, float, str]] = {
    "DOR (4NCG)": ("4NCG", "2KW", "X-ray", 2.58, "no"),
    "DOR (7Z2G)": ("7Z2G", "2KW", "cryo-EM", 3.65, "38-mer DNA"),
    "RPV (2ZD1)": ("2ZD1", "T27", "X-ray", 1.80, "no"),
    "EFV (1FK9)": ("1FK9", "EFZ", "X-ray", 2.50, "no"),
    "NVP (1VRT)": ("1VRT", "NVP", "X-ray", 2.20, "no"),
}
REFERENCE = "DOR (4NCG)"

# Tyr chi1 is N-CA-CB-CG, chi2 is CA-CB-CG-CD1.
CHI1 = ("N", "CA", "CB", "CG")
CHI2 = ("CA", "CB", "CG", "CD1")
RING = ("CG", "CD1", "CD2", "CE1", "CE2", "CZ")
TYROSINES = {"Tyr181": 181, "Tyr188": 188}


def chi(res: dict[str, np.ndarray], names: tuple[str, ...]) -> float:
    if not all(n in res for n in names):
        return float("nan")
    return dihedral(np.array([res[n] for n in names]))


def ring_centroid(res: dict[str, np.ndarray]) -> np.ndarray:
    present = [res[n] for n in RING if n in res]
    if len(present) < 4:
        return np.full(3, np.nan)
    return np.mean(present, axis=0)


def min_distance(a: dict[str, np.ndarray], b: dict[str, np.ndarray]) -> float:
    if not a or not b:
        return float("nan")
    A = np.array(list(a.values()))
    B = np.array(list(b.values()))
    return float(np.linalg.norm(A[:, None, :] - B[None, :, :], axis=-1).min())


def circular_delta(a: float, b: float) -> float:
    """Signed smallest difference a - b, in degrees, on (-180, 180]."""
    if np.isnan(a) or np.isnan(b):
        return float("nan")
    d = (a - b + 180.0) % 360.0 - 180.0
    return float(d)


def main(argv: list[str] | None = None) -> int:
    root = Path(__file__).resolve().parents[3]
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--structure-dir", type=Path,
                    default=root / "data/structures/apo_reference")
    ap.add_argument("--output-dir", type=Path,
                    default=root / "results/analysis/nnrti_pose_comparison")
    args = ap.parse_args(argv)

    tables = args.output_dir / "tables"
    tables.mkdir(parents=True, exist_ok=True)

    loaded: dict[str, Struct] = {}
    for label, (pdb_id, ligand, *_rest) in STRUCTURES.items():
        path = args.structure_dir / f"{pdb_id}.cif"
        if not path.exists():
            print(f"[skip] {label}: {path.name} not found")
            continue
        loaded[label] = Struct(path, label, ligand=ligand)

    if REFERENCE not in loaded:
        raise SystemExit(f"reference {REFERENCE} missing")
    ref = loaded[REFERENCE]

    # the pocket frame: NNIBP-lining Calpha atoms present in every structure
    want = [("A", r) for r in NNIBP_P66] + [("B", r) for r in NNIBP_P51]
    keys = [k for k in want if all(k in s.res and "CA" in s.res[k] for s in loaded.values())]
    print(f"pocket frame: {len(keys)} of {len(want)} NNIBP Calpha atoms common to all structures")

    rows = []
    per_atom = []
    for label, st in loaded.items():
        pdb_id, ligand, method, resolution, na = STRUCTURES[label]

        # superpose this structure's pocket onto the reference pocket
        rt = kabsch(st.ca(keys), ref.ca(keys))
        fit_rmsd = rmsd(apply(rt, st.ca(keys)), ref.ca(keys))

        row = {
            "structure": label, "pdb_id": pdb_id, "ligand": ligand,
            "method": method, "resolution_A": resolution, "nucleic_acid": na,
            "pocket_ca_rmsd_vs_ref_A": round(fit_rmsd, 3),
        }

        for name, resid in TYROSINES.items():
            key = ("A", resid)
            res = st.res.get(key, {})
            refres = ref.res.get(key, {})
            mutated = st.resname.get(key, "???") != "TYR"
            row[f"{name}_residue"] = st.resname.get(key, "???")
            if mutated or not res:
                for suffix in ("chi1_deg", "chi2_deg", "dchi1_vs_ref_deg",
                               "ring_shift_vs_ref_A", "min_dist_to_ligand_A"):
                    row[f"{name}_{suffix}"] = float("nan")
                continue

            c1, c2 = chi(res, CHI1), chi(res, CHI2)
            row[f"{name}_chi1_deg"] = round(c1, 1)
            row[f"{name}_chi2_deg"] = round(c2, 1)
            row[f"{name}_dchi1_vs_ref_deg"] = round(circular_delta(c1, chi(refres, CHI1)), 1)

            # ring centroid displacement, in the common pocket frame
            moved = {n: apply(rt, v[None, :])[0] for n, v in res.items()}
            shift = np.linalg.norm(ring_centroid(moved) - ring_centroid(refres))
            row[f"{name}_ring_shift_vs_ref_A"] = round(float(shift), 2)

            # contact with that structure's OWN ligand -- no superposition
            row[f"{name}_min_dist_to_ligand_A"] = round(min_distance(res, st.lig), 2)

            for atom, xyz in sorted(moved.items()):
                if atom in RING or atom in ("N", "CA", "CB", "OH"):
                    per_atom.append({"structure": label, "residue": name, "atom": atom,
                                     "x": round(float(xyz[0]), 3),
                                     "y": round(float(xyz[1]), 3),
                                     "z": round(float(xyz[2]), 3)})

        # write the superposed coordinates so a render and the numbers come
        # from one and the same fit
        struct_dir = args.output_dir / "structures"
        struct_dir.mkdir(parents=True, exist_ok=True)
        import gemmi
        gst = gemmi.read_structure(str(args.structure_dir / f"{pdb_id}.cif"))
        gst.setup_entities(); gst.remove_alternative_conformations()
        gst.remove_hydrogens(); gst.remove_waters()
        R, t = rt
        for model in gst:
            for chain in model:
                for res in chain:
                    for atom in res:
                        v = np.array([atom.pos.x, atom.pos.y, atom.pos.z])
                        w = R @ v + t
                        atom.pos = gemmi.Position(*w)
        gst.setup_entities()
        gst.write_pdb(str(struct_dir / f"{pdb_id}_on_4NCG_nnibp_frame.pdb"))

        pdb = ensure_pdb(args.structure_dir, pdb_id)
        try:
            row["nnibp_volume_A3"] = round(pocket_volume(pdb), 1)
        except Exception as exc:
            row["nnibp_volume_A3"] = float("nan")
            print(f"[warn] {label}: pocket volume failed ({exc})")

        rows.append(row)

    df = pd.DataFrame(rows)
    out = tables / "nnrti_tyrosine_comparison.csv"
    df.to_csv(out, index=False)
    pd.DataFrame(per_atom).to_csv(tables / "nnrti_tyrosine_atoms.csv", index=False)

    show = ["structure", "resolution_A",
            "Tyr181_chi1_deg", "Tyr181_ring_shift_vs_ref_A", "Tyr181_min_dist_to_ligand_A",
            "Tyr188_chi1_deg", "Tyr188_ring_shift_vs_ref_A", "Tyr188_min_dist_to_ligand_A",
            "nnibp_volume_A3"]
    print()
    print(df[show].to_string(index=False))
    print(f"\nwrote {out}")

    summary = {
        "reference": REFERENCE,
        "pocket_frame_atoms": len(keys),
        "structures": {r["structure"]: {k: v for k, v in r.items() if k != "structure"}
                       for r in rows},
    }
    (tables / "nnrti_tyrosine_comparison.json").write_text(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
