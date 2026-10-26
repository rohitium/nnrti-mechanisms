#!/usr/bin/env python3
"""Compare the two independent experimental RT-DOR structures with each other.

Reviewer 1 questions whether a cycle built on one DOR-bound crystal structure
can stand for the biologically relevant RT-DOR ensemble. One reading of that
objection is specifically about 4NCG: that the pose we sample around is a
crystallographic artefact rather than the pose DOR actually adopts.

That reading is testable. A second RT-DOR structure exists, determined by a
different method, in a different laboratory, with nucleic acid bound:

    4NCG   2.58 A X-ray, RT + DOR, no nucleic acid        (Lai 2014)
    7Z2G   3.65 A cryo-EM, RT + DOR + 38-mer DNA aptaper  (2022)

If the two agree on where DOR sits, which NNIBP side chains line it, and how
they are oriented, then the 4NCG pose is not an artefact of one crystal form,
and it survives the addition of nucleic acid -- which is the physiological
state Reviewer 2 asks about separately.

Both ligands are the same chemical component (2KW), so ligand atoms are matched
by name. Residues are matched on canonical HIV-1 RT numbering, which both
depositions use for p66 (chain A) and p51 (chain B); the -3 offset of the
simulation topology is applied only where this module calls into
nnrti.analysis.metrics, which expects topology numbering.

Three superpositions are reported, because the answer depends on the frame:

    whole    both subunits' Ca -- whether the enzymes are the same shape at all
    palm     p66 1-315 Ca      -- a frame that excludes the thumb, which moves
                                  when DNA binds and is not what is at issue
    nnibp    the 16 NNIBP Ca   -- the pocket's own frame, in which the ligand
                                  RMSD answers the reviewer's question directly

The contact and gate distances at the end are superposition-free: they do not
depend on any fit, so they cannot be flattered by the choice of reference set.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.compare_dor_crystal_poses
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

AA = set(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR "
    "TRP TYR VAL".split()
)

#: canonical HIV-1 RT numbering, Cilento, Kirby & Sarafianos
NNIBP_P66 = (100, 101, 103, 106, 107, 108, 179, 181, 188, 189, 190, 227, 229, 234, 318)
NNIBP_P51 = (138,)

#: p66 subdomains, canonical numbering (Kohlstaedt 1992 boundaries). Reported
#: because the two depositions must differ somewhere -- one has DNA bound and
#: the other does not -- and the question is whether they differ at the NNIBP
#: or at the domains that grip nucleic acid.
SUBDOMAINS = (
    ("p66 fingers", "A", ((1, 85), (118, 155))),
    ("p66 palm", "A", ((86, 117), (156, 236))),
    ("p66 thumb", "A", ((237, 318),)),
    ("p66 connection", "A", ((319, 426),)),
    ("p66 RNase H", "A", ((427, 560),)),
    ("p51", "B", ((1, 440),)),
)

#: p66 polymerase subdomains, canonical numbering
PALM_RANGE = (1, 315)

#: the four NNIBP gates measured in the tunnel-dynamics analysis
GATES = (
    ("K101-Y188", ("A", 101), ("A", 188)),
    ("V106-Y181", ("A", 106), ("A", 181)),
    ("L100-F227", ("A", 100), ("A", 227)),
    ("K103-P225", ("A", 103), ("A", 225)),
)

#: DOR rotatable dihedrals in chemical-component (2KW) atom names. Derived by
#: graph isomorphism against the simulation topology names used in
#: nnrti.cli.compute_dor_torsions, so the values are directly comparable with
#: dor_torsions_per_frame.csv.
TORSIONS = (
    ("tau1", ("C12", "C7", "O", "C2")),
    ("tau2", ("C7", "O", "C2", "C1")),
    ("tau3", ("C10", "N11", "C17", "C18")),
    ("tau4", ("N11", "C17", "C18", "N22")),
)

CHI1 = {
    "ARG": ("N", "CA", "CB", "CG"), "ASN": ("N", "CA", "CB", "CG"),
    "ASP": ("N", "CA", "CB", "CG"), "CYS": ("N", "CA", "CB", "SG"),
    "GLN": ("N", "CA", "CB", "CG"), "GLU": ("N", "CA", "CB", "CG"),
    "HIS": ("N", "CA", "CB", "CG"), "ILE": ("N", "CA", "CB", "CG1"),
    "LEU": ("N", "CA", "CB", "CG"), "LYS": ("N", "CA", "CB", "CG"),
    "MET": ("N", "CA", "CB", "CG"), "PHE": ("N", "CA", "CB", "CG"),
    "PRO": ("N", "CA", "CB", "CG"), "SER": ("N", "CA", "CB", "OG"),
    "THR": ("N", "CA", "CB", "OG1"), "TRP": ("N", "CA", "CB", "CG"),
    "TYR": ("N", "CA", "CB", "CG"), "VAL": ("N", "CA", "CB", "CG1"),
}
CHI2 = {
    "ARG": ("CA", "CB", "CG", "CD"), "ASN": ("CA", "CB", "CG", "OD1"),
    "ASP": ("CA", "CB", "CG", "OD1"), "GLN": ("CA", "CB", "CG", "CD"),
    "GLU": ("CA", "CB", "CG", "CD"), "HIS": ("CA", "CB", "CG", "ND1"),
    "ILE": ("CA", "CB", "CG1", "CD1"), "LEU": ("CA", "CB", "CG", "CD1"),
    "LYS": ("CA", "CB", "CG", "CD"), "MET": ("CA", "CB", "CG", "SD"),
    "PHE": ("CA", "CB", "CG", "CD1"), "TRP": ("CA", "CB", "CG", "CD1"),
    "TYR": ("CA", "CB", "CG", "CD1"),
}

BACKBONE = {"N", "CA", "C", "O", "OXT"}

#: side chains whose terminal atoms are indistinguishable. A ring flip or a
#: carboxylate swap is the same structure, so an atom-name-matched RMSD has to
#: try both assignments and keep the better one, or Phe227 appears to move 2 A
#: when its ring has merely been modelled the other way round.
SYMMETRIC = {
    "PHE": (("CD1", "CD2"), ("CE1", "CE2")),
    "TYR": (("CD1", "CD2"), ("CE1", "CE2")),
    "ASP": (("OD1", "OD2"),),
    "GLU": (("OE1", "OE2"),),
    "ARG": (("NH1", "NH2"),),
}


# --------------------------------------------------------------------------- #
# geometry
# --------------------------------------------------------------------------- #
def kabsch(mobile: np.ndarray, target: np.ndarray):
    """Return (rotation, translation) that best maps `mobile` onto `target`."""
    mc, tc = mobile.mean(0), target.mean(0)
    h = (mobile - mc).T @ (target - tc)
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return rot, tc - rot @ mc


def apply(rt, xyz: np.ndarray) -> np.ndarray:
    rot, trans = rt
    return (rot @ xyz.T).T + trans


def rmsd(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(((a - b) ** 2).sum(1).mean()))


def dihedral(p: np.ndarray) -> float:
    b0, b1, b2 = p[0] - p[1], p[2] - p[1], p[3] - p[2]
    b1 /= np.linalg.norm(b1)
    v = b0 - (b0 @ b1) * b1
    w = b2 - (b2 @ b1) * b1
    return float(np.degrees(np.arctan2(np.cross(b1, v) @ w, v @ w)))


def circ_diff(a: float, b: float) -> float:
    """Signed smallest difference a - b, in degrees, on the circle."""
    return float((a - b + 180.0) % 360.0 - 180.0)


def circ_diff_180(a: float, b: float) -> float:
    """As circ_diff, but on a two-fold symmetric angle: wrapped to +/-90."""
    return float((a - b + 90.0) % 180.0 - 90.0)


def symmetry_corrected_rmsd(
    ref_atoms: dict[str, np.ndarray],
    mob_atoms: dict[str, np.ndarray],
    names: list[str],
    resname: str,
) -> tuple[float, bool]:
    """RMSD over `names`, minimised over the side chain's equivalent atoms."""
    pairs = [p for p in SYMMETRIC.get(resname, ())
             if p[0] in names and p[1] in names]
    a = np.array([ref_atoms[n] for n in names])
    plain = rmsd(a, np.array([mob_atoms[n] for n in names]))
    if not pairs:
        return plain, False
    swapped = dict(mob_atoms)
    for x, y in pairs:
        swapped[x], swapped[y] = mob_atoms[y], mob_atoms[x]
    flipped = rmsd(a, np.array([swapped[n] for n in names]))
    return (min(plain, flipped), flipped < plain)


def _ligand_graph(atoms: dict[str, np.ndarray]):
    import networkx as nx

    cov = {"C": 0.77, "N": 0.75, "O": 0.73, "F": 0.71, "CL": 0.99, "S": 1.02}
    names = list(atoms)
    el = {n: ("CL" if n.upper().startswith("CL") else n[0].upper()) for n in names}
    g = nx.Graph()
    for n in names:
        g.add_node(n, el=el[n])
    for i, n in enumerate(names):
        for m in names[i + 1:]:
            if float(np.linalg.norm(atoms[n] - atoms[m])) < (
                    cov[el[n]] + cov[el[m]]) * 1.3:
                g.add_edge(n, m)
    return g


def ligand_automorphisms(atoms: dict[str, np.ndarray]) -> list[dict[str, str]]:
    """Atom-name permutations that leave the ligand's bond graph unchanged.

    DOR's trifluoromethyl fluorines are indistinguishable, so an atom-name
    matched RMSD charges the pyridinone up to half an angstrom for a methyl
    rotation that is not a conformational difference at all. Minimising over
    the automorphism group removes that.
    """
    from networkx.algorithms.isomorphism import GraphMatcher

    g = _ligand_graph(atoms)
    gm = GraphMatcher(g, g, node_match=lambda a, b: a["el"] == b["el"])
    return [dict(m) for m in gm.isomorphisms_iter()]


def ring_systems(atoms: dict[str, np.ndarray]) -> dict[str, list[str]]:
    """Partition DOR heavy atoms into its three rings plus the linker.

    Rings come from bond connectivity, as in compute_dor_moiety_contacts, and
    are named by composition rather than by atom name: the five-membered ring
    with three nitrogens is the triazolinone, the six-ring with one nitrogen is
    the pyridinone, the all-carbon six-ring is the chlorocyanophenyl. The ether
    oxygen and the methylene carbon are the linker; every other substituent
    joins the ring it is bonded to.
    """
    import networkx as nx

    cov = {"C": 0.77, "N": 0.75, "O": 0.73, "F": 0.71, "CL": 0.99, "S": 1.02}
    names = list(atoms)
    el = {n: ("CL" if n.upper().startswith("CL") else n[0].upper()) for n in names}
    g = nx.Graph()
    g.add_nodes_from(names)
    for i, n in enumerate(names):
        for m in names[i + 1:]:
            d = float(np.linalg.norm(atoms[n] - atoms[m]))
            if d < (cov[el[n]] + cov[el[m]]) * 1.3:
                g.add_edge(n, m)

    rings: list[set[str]] = []
    for cycle in nx.cycle_basis(g):
        rings.append(set(cycle))
    out: dict[str, list[str]] = {}
    for ring in rings:
        n_n = sum(1 for a in ring if el[a] == "N")
        if len(ring) == 5 and n_n == 3:
            key = "triazolinone"
        elif len(ring) == 6 and n_n == 1:
            key = "pyridinone"
        elif len(ring) == 6 and n_n == 0:
            key = "chlorocyanophenyl"
        else:
            key = f"ring_{len(ring)}"
        out[key] = sorted(ring)

    in_ring = {a for ring in rings for a in ring}
    linker, assigned = [], {}
    for a in names:
        if a in in_ring:
            continue
        nbrs = [b for b in g.neighbors(a)]
        ring_nbrs = [b for b in nbrs if b in in_ring]
        if el[a] == "O" and len(ring_nbrs) == 2:
            linker.append(a)                      # ether oxygen
        elif el[a] == "C" and len(ring_nbrs) == 1 and len(nbrs) == 2:
            linker.append(a)                      # methylene
        else:
            assigned[a] = ring_nbrs
    for a, ring_nbrs in assigned.items():
        home = None
        for key, members in out.items():
            if any(b in members for b in ring_nbrs):
                home = key
                break
        if home is None:                          # two bonds out from a ring
            for key, members in out.items():
                if any(b in members
                       for nb in g.neighbors(a) for b in g.neighbors(nb)):
                    home = key
                    break
        if home is not None:
            out[home] = sorted(out[home] + [a])
        else:
            linker.append(a)
    out["linker"] = sorted(linker)
    return out


# --------------------------------------------------------------------------- #
# structure access
# --------------------------------------------------------------------------- #
class Struct:
    """Residue- and atom-level access to one deposition, by canonical numbering."""

    def __init__(self, path: Path, label: str, ligand: str = "2KW"):
        import gemmi

        st = gemmi.read_structure(str(path))
        st.setup_entities()
        st.remove_alternative_conformations()
        st.remove_hydrogens()
        st.remove_waters()
        self.label = label
        self.path = path
        self.res: dict[tuple[str, int], dict[str, np.ndarray]] = {}
        self.resname: dict[tuple[str, int], str] = {}
        self.lig: dict[str, np.ndarray] = {}
        for chain in st[0]:
            for r in chain:
                atoms = {
                    a.name: np.array([a.pos.x, a.pos.y, a.pos.z])
                    for a in r
                    if a.element.name != "H"
                }
                if r.name == ligand and chain.name == "A":
                    self.lig = atoms
                elif r.name in AA:
                    key = (chain.name, r.seqid.num)
                    self.res[key] = atoms
                    self.resname[key] = r.name
        if not self.lig:
            raise SystemExit(f"{label}: ligand {ligand} not found in chain A")

    def ca(self, keys) -> np.ndarray:
        return np.array([self.res[k]["CA"] for k in keys])

    def chi(self, key, table) -> float:
        names = table.get(self.resname[key])
        if not names:
            return float("nan")
        atoms = self.res[key]
        if any(n not in atoms for n in names):
            return float("nan")
        return dihedral(np.array([atoms[n] for n in names]))

    def heavy(self, key, sidechain_only: bool = False) -> dict[str, np.ndarray]:
        atoms = self.res[key]
        if not sidechain_only:
            return atoms
        return {n: v for n, v in atoms.items() if n not in BACKBONE}

    def ligand_xyz(self, names) -> np.ndarray:
        return np.array([self.lig[n] for n in names])

    def min_dist_to_ligand(self, key) -> float:
        lig = np.array(list(self.lig.values()))
        res = np.array(list(self.res[key].values()))
        return float(np.linalg.norm(res[:, None, :] - lig[None, :, :], axis=2).min())

    def min_dist(self, key_a, key_b) -> float:
        a = np.array(list(self.res[key_a].values()))
        b = np.array(list(self.res[key_b].values()))
        return float(np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2).min())


# --------------------------------------------------------------------------- #
def ensure_pdb(structure_dir: Path, pdb_id: str) -> Path:
    """Return the PDB file for `pdb_id`, writing it from the mmCIF if absent.

    Only the mmCIF depositions are kept in the repository -- .pdb is ignored --
    so a fresh checkout has the CIF and not the PDB. MDAnalysis reads the PDB
    more reliably than the CIF here, so convert rather than depend on one.
    """
    import gemmi

    pdb = structure_dir / f"{pdb_id}.pdb"
    if pdb.exists():
        return pdb
    cif = structure_dir / f"{pdb_id}.cif"
    if not cif.exists():
        raise SystemExit(f"neither {pdb.name} nor {cif.name} in {structure_dir}")
    st = gemmi.read_structure(str(cif))
    st.setup_entities()
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    st.write_pdb(str(pdb))
    return pdb


def pocket_volume(path: Path, offset: int = -3) -> float:
    """NNIBP pocket volume through the same routine used on the trajectories."""
    import MDAnalysis as mda

    from nnrti.analysis.metrics import pocket_volume_proxy_from_universe

    u = mda.Universe(str(path))
    u.residues.resids = u.residues.resids + offset
    return float(pocket_volume_proxy_from_universe(u, ""))


def main(argv: list[str] | None = None) -> int:
    root = Path(__file__).resolve().parents[3]
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--structure-dir", type=Path,
                    default=root / "data/structures/apo_reference")
    ap.add_argument("--reference", default="4NCG")
    ap.add_argument("--mobile", default="7Z2G")
    ap.add_argument("--torsion-frames", type=Path,
                    default=root / "results/analysis/dor_torsions/tables/"
                                   "dor_torsions_per_frame.csv")
    ap.add_argument("--output-dir", type=Path,
                    default=root / "results/analysis/crystal_pose_comparison")
    args = ap.parse_args(argv)

    tables = args.output_dir / "tables"
    tables.mkdir(parents=True, exist_ok=True)

    ref = Struct(ensure_pdb(args.structure_dir, args.reference), args.reference)
    mob = Struct(ensure_pdb(args.structure_dir, args.mobile), args.mobile)

    # ---- residue sets both depositions resolve ---------------------------- #
    common = sorted(
        k for k in ref.res
        if k in mob.res
        and ref.resname[k] == mob.resname[k]
        and "CA" in ref.res[k] and "CA" in mob.res[k]
    )
    p66 = [k for k in common if k[0] == "A"]
    p51 = [k for k in common if k[0] == "B"]
    palm = [k for k in p66 if PALM_RANGE[0] <= k[1] <= PALM_RANGE[1]]
    nnibp = [("A", r) for r in NNIBP_P66] + [("B", r) for r in NNIBP_P51]
    missing = [k for k in nnibp if k not in common]
    if missing:
        raise SystemExit(f"NNIBP residues unresolved or mismatched: {missing}")

    # ---- three superpositions --------------------------------------------- #
    frames = {}
    fit_rows = []
    for name, keys in (("whole", common), ("palm", palm), ("nnibp", nnibp)):
        rt = kabsch(mob.ca(keys), ref.ca(keys))
        frames[name] = rt
        fit_rows.append({
            "frame": name,
            "n_fit_residues": len(keys),
            "fit_ca_rmsd_A": round(rmsd(apply(rt, mob.ca(keys)), ref.ca(keys)), 3),
            "p66_ca_rmsd_A": round(rmsd(apply(rt, mob.ca(p66)), ref.ca(p66)), 3),
            "p51_ca_rmsd_A": round(rmsd(apply(rt, mob.ca(p51)), ref.ca(p51)), 3),
            "nnibp_ca_rmsd_A": round(rmsd(apply(rt, mob.ca(nnibp)), ref.ca(nnibp)), 3),
        })
    fit_df = pd.DataFrame(fit_rows)

    # where do the two structures differ? fit the palm, which neither DNA
    # binding nor NNRTI binding is expected to move, and look at each subdomain
    rt_palm = frames["palm"]
    sub_rows = []
    for name, chain, spans in SUBDOMAINS:
        keys = [k for k in common
                if k[0] == chain and any(lo <= k[1] <= hi for lo, hi in spans)]
        if not keys:
            continue
        sub_rows.append({
            "subdomain": name,
            "n_residues": len(keys),
            "ca_rmsd_palm_frame_A": round(
                rmsd(apply(rt_palm, mob.ca(keys)), ref.ca(keys)), 3),
        })
    sub_rows.append({
        "subdomain": "p66 NNIBP lining (15 residues)",
        "n_residues": len(NNIBP_P66),
        "ca_rmsd_palm_frame_A": round(rmsd(
            apply(rt_palm, mob.ca([("A", r) for r in NNIBP_P66])),
            ref.ca([("A", r) for r in NNIBP_P66])), 3),
    })
    sub_df = pd.DataFrame(sub_rows)

    # ---- ligand ----------------------------------------------------------- #
    lig_names = sorted(set(ref.lig) & set(mob.lig))
    if len(lig_names) != len(ref.lig) or len(lig_names) != len(mob.lig):
        raise SystemExit("ligand atom names differ between depositions")
    # resolve the trifluoromethyl's indistinguishable fluorines in the pocket
    # frame, then keep that labelling for every ligand number below
    autos = ligand_automorphisms(ref.lig)
    ref_lig = ref.ligand_xyz(lig_names)
    best = min(
        autos,
        key=lambda m: rmsd(
            ref_lig,
            apply(frames["nnibp"], np.array([mob.lig[m[n]] for n in lig_names])),
        ),
    )
    n_relabelled = sum(1 for n in lig_names if best[n] != n)
    mob.lig = {n: mob.lig[best[n]] for n in lig_names}
    mob_lig = mob.ligand_xyz(lig_names)

    lig_rows = []
    for name, rt in frames.items():
        moved = apply(rt, mob_lig)
        lig_rows.append({
            "frame": name,
            "n_ligand_atoms": len(lig_names),
            "ligand_rmsd_A": round(rmsd(moved, ref_lig), 3),
            "ligand_centroid_shift_A": round(
                float(np.linalg.norm(moved.mean(0) - ref_lig.mean(0))), 3),
        })
    rt_lig = kabsch(mob_lig, ref_lig)
    lig_rows.append({
        "frame": "ligand-on-ligand (internal conformation only)",
        "n_ligand_atoms": len(lig_names),
        "ligand_rmsd_A": round(rmsd(apply(rt_lig, mob_lig), ref_lig), 3),
        "ligand_centroid_shift_A": 0.0,
    })
    lig_df = pd.DataFrame(lig_rows)

    # per-moiety, in the pocket's frame: a whole-ligand RMSD cannot say whether
    # one ring system is placed differently from the rest
    moieties = ring_systems(ref.lig)
    moi_rows = []
    for name in ("chlorocyanophenyl", "pyridinone", "triazolinone", "linker"):
        names = [n for n in moieties.get(name, []) if n in mob.lig]
        if not names:
            continue
        a = ref.ligand_xyz(names)
        b = apply(frames["nnibp"], mob.ligand_xyz(names))
        moi_rows.append({
            "moiety": name,
            "n_atoms": len(names),
            "atoms": " ".join(names),
            "rmsd_nnibp_frame_A": round(rmsd(a, b), 3),
            "centroid_shift_A": round(
                float(np.linalg.norm(b.mean(0) - a.mean(0))), 3),
        })
    moi_df = pd.DataFrame(moi_rows)

    # ---- per-residue, in the NNIBP frame ---------------------------------- #
    rt_n = frames["nnibp"]
    res_rows = []
    for key in nnibp:
        chain, resid = key
        ref_at, mob_at = ref.heavy(key), mob.heavy(key)
        resname = ref.resname[key]
        shared = sorted(set(ref_at) & set(mob_at))
        sc = [n for n in shared if n not in BACKBONE]
        moved = {n: v for n, v in zip(shared, apply(
            rt_n, np.array([mob_at[n] for n in shared])))}
        all_rmsd, _ = symmetry_corrected_rmsd(ref_at, moved, shared, resname)
        sc_rmsd, flipped = (symmetry_corrected_rmsd(ref_at, moved, sc, resname)
                            if sc else (np.nan, False))
        two_fold = resname in SYMMETRIC
        dchi2 = (circ_diff_180 if two_fold else circ_diff)(
            mob.chi(key, CHI2), ref.chi(key, CHI2))
        res_rows.append({
            "subunit": "p66" if chain == "A" else "p51",
            "resid": resid,
            "resname": resname,
            "ca_deviation_A": round(float(np.linalg.norm(
                apply(rt_n, mob_at["CA"]) - ref_at["CA"])), 3),
            "all_heavy_rmsd_A": round(all_rmsd, 3),
            "sidechain_rmsd_A": round(sc_rmsd, 3) if sc else np.nan,
            "terminal_atoms_swapped": flipped,
            "n_atoms_compared": len(shared),
            f"chi1_{args.reference}_deg": round(ref.chi(key, CHI1), 1),
            f"chi1_{args.mobile}_deg": round(mob.chi(key, CHI1), 1),
            "delta_chi1_deg": round(circ_diff(mob.chi(key, CHI1),
                                              ref.chi(key, CHI1)), 1),
            "delta_chi2_deg": round(dchi2, 1),
            "chi2_two_fold": two_fold,
            f"min_dist_to_DOR_{args.reference}_A": round(
                ref.min_dist_to_ligand(key), 2),
            f"min_dist_to_DOR_{args.mobile}_A": round(
                mob.min_dist_to_ligand(key), 2),
        })
    res_df = pd.DataFrame(res_rows)
    res_df["delta_min_dist_A"] = (
        res_df[f"min_dist_to_DOR_{args.mobile}_A"]
        - res_df[f"min_dist_to_DOR_{args.reference}_A"]
    ).round(2)

    # ---- gates and torsions (superposition-free) -------------------------- #
    gate_rows = []
    for name, a, b in GATES:
        if a not in ref.res or b not in ref.res or a not in mob.res or b not in mob.res:
            continue
        gate_rows.append({
            "gate": name,
            f"{args.reference}_min_heavy_dist_A": round(ref.min_dist(a, b), 2),
            f"{args.mobile}_min_heavy_dist_A": round(mob.min_dist(a, b), 2),
        })
    gate_df = pd.DataFrame(gate_rows)
    gate_df["delta_A"] = (gate_df[f"{args.mobile}_min_heavy_dist_A"]
                          - gate_df[f"{args.reference}_min_heavy_dist_A"]).round(2)

    tor_rows = []
    wt = None
    if args.torsion_frames.exists():
        per_frame = pd.read_csv(args.torsion_frames)
        wt = per_frame[per_frame["safe_label"].str.lower() == "wt"]
        if wt.empty:
            wt = None
    for name, atoms in TORSIONS:
        row = {"torsion": name, "atoms": "-".join(atoms)}
        for st in (ref, mob):
            row[f"{st.label}_deg"] = round(
                dihedral(np.array([st.lig[a] for a in atoms])), 1)
        row["delta_deg"] = round(
            circ_diff(row[f"{mob.label}_deg"], row[f"{ref.label}_deg"]), 1)
        if wt is not None:
            vals = wt[name].to_numpy()
            r = np.deg2rad(vals)
            mean = float(np.rad2deg(np.arctan2(np.sin(r).mean(), np.cos(r).mean())))
            resultant = float(np.hypot(np.sin(r).mean(), np.cos(r).mean()))
            sd = float(np.rad2deg(np.sqrt(-2.0 * np.log(max(resultant, 1e-12)))))
            row["WT_MD_circmean_deg"] = round(mean, 1)
            row["WT_MD_circsd_deg"] = round(sd, 1)
            for st in (ref, mob):
                row[f"{st.label}_z_in_WT_MD"] = round(
                    circ_diff(row[f"{st.label}_deg"], mean) / sd, 2) if sd else np.nan
        tor_rows.append(row)
    tor_df = pd.DataFrame(tor_rows)

    # ---- pocket volume ---------------------------------------------------- #
    vol_rows = []
    for st, method in ((ref, "X-ray 2.58 A"), (mob, "cryo-EM 3.65 A")):
        vol_rows.append({
            "structure": st.label,
            "method": method,
            "nucleic_acid": "no" if st.label == "4NCG" else "38-mer DNA",
            "pocket_volume_A3": round(pocket_volume(st.path), 1),
        })
    vol_df = pd.DataFrame(vol_rows)

    # ---- the superposed coordinates, so a figure uses this exact fit ------ #
    import gemmi

    struct_dir = args.output_dir / "structures"
    struct_dir.mkdir(parents=True, exist_ok=True)
    rot, trans = frames["nnibp"]
    st = gemmi.read_structure(str(mob.path))
    st.setup_entities()
    for model in st:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    v = rot @ np.array([atom.pos.x, atom.pos.y, atom.pos.z]) + trans
                    atom.pos = gemmi.Position(*v)
    out_pdb = struct_dir / f"{args.mobile}_on_{args.reference}_nnibp_frame.pdb"
    st.write_pdb(str(out_pdb))

    # ---- unresolved residues near the pocket (a caveat, not a result) ----- #
    centre = ref.ca(nnibp).mean(0)
    near_ref = {k for k in ref.res
                if np.linalg.norm(ref.res[k]["CA"] - centre) <= 15.0}
    gaps = sorted(k for k in near_ref if k not in mob.res)
    truncated = sorted(
        k for k in near_ref & set(mob.res)
        if set(ref.heavy(k, True)) - set(mob.heavy(k, True))
    )

    for name, df in (("superposition_rmsd", fit_df), ("ligand_rmsd", lig_df),
                     ("ligand_moiety_rmsd", moi_df),
                     ("subdomain_rmsd", sub_df),
                     ("nnibp_residue_comparison", res_df), ("nnibp_gates", gate_df),
                     ("dor_torsions_crystal", tor_df), ("pocket_volume", vol_df)):
        df.to_csv(tables / f"{name}.csv", index=False)

    d_lig = res_df["delta_min_dist_A"]
    nn = lig_df.set_index("frame")
    summary = f"""# {args.reference} vs {args.mobile}: is the DOR pose structure-specific?

{args.reference}: X-ray 2.58 A, RT + DOR, no nucleic acid.
{args.mobile}: cryo-EM 3.65 A, RT + DOR + 38-mer DNA aptamer, independent group.
Ligand is chemical component 2KW in both; {len(lig_names)} heavy atoms matched by
name, up to the {len(autos)} automorphisms of its bond graph (the CF3 fluorines);
{n_relabelled} atom labels were permuted to resolve those.
{len(common)} residues resolved in both ({len(p66)} p66, {len(p51)} p51).

## Enzyme

{fit_df.to_markdown(index=False)}

Where the two differ, with the palm as the reference frame:

{sub_df.to_markdown(index=False)}

## DOR

{lig_df.to_markdown(index=False)}

Ligand heavy-atom RMSD in the pocket's own frame: \
{nn.loc['nnibp', 'ligand_rmsd_A']:.2f} A over {len(lig_names)} atoms.
Internal conformation alone (ligand fitted on ligand): \
{nn.loc['ligand-on-ligand (internal conformation only)', 'ligand_rmsd_A']:.2f} A.

### By moiety, in the NNIBP frame

{moi_df.drop(columns=['atoms']).to_markdown(index=False)}

### Rotatable dihedrals

{tor_df.to_markdown(index=False)}

## NNIBP residues, in the NNIBP frame

{res_df.to_markdown(index=False)}

Side-chain RMSD: median {res_df['sidechain_rmsd_A'].median():.2f} A, \
max {res_df['sidechain_rmsd_A'].max():.2f} A \
({res_df.loc[res_df['sidechain_rmsd_A'].idxmax(), 'resname']}\
{res_df.loc[res_df['sidechain_rmsd_A'].idxmax(), 'resid']}).
|delta chi1| > 30 deg: {int((res_df['delta_chi1_deg'].abs() > 30).sum())} of \
{int(res_df['delta_chi1_deg'].notna().sum())} residues.

## Superposition-free checks

Ligand-to-residue minimum heavy-atom distances agree to a median of \
{d_lig.abs().median():.2f} A (max {d_lig.abs().max():.2f} A); Pearson r = \
{res_df[f'min_dist_to_DOR_{args.reference}_A'].corr(res_df[f'min_dist_to_DOR_{args.mobile}_A']):.3f}.

{gate_df.to_markdown(index=False)}

{vol_df.to_markdown(index=False)}

## Caveats

{args.mobile} is 3.65 A, so its side-chain positions carry more uncertainty than
{args.reference}'s. Within 15 A of the NNIBP Ca centroid, {args.mobile} leaves
{len(gaps)} residues unmodelled{' (' + ', '.join(f'{ref.resname[k]}{k[1]}' for k in gaps) + ')' if gaps else ''}
and truncates the side chain of {len(truncated)} more\
{' (' + ', '.join(f'{ref.resname[k]}{k[1]}' for k in truncated) + ')' if truncated else ''}.
Those residues are excluded atom-by-atom from the per-residue comparison; the
pocket-volume figures are not corrected for them, so {args.mobile}'s volume is
an upper bound.
"""
    (args.output_dir / "summary.md").write_text(summary)
    print(summary)
    print(f"written: {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
