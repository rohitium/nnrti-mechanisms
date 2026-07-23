"""RDKit/OpenFF-backed stand-ins for the OpenEye APIs Perses expects.

Perses point-mutation setup imports ``openeye.oechem`` unconditionally and uses
OpenEye for ligand SDF parsing, amino-acid template graphs, and MCS atom mapping.
This module installs lightweight replacements so hybrid prep can run without a
licensed OpenEye install.
"""

from __future__ import annotations

import itertools
import sys
import types
from pathlib import Path

import numpy as np
from openmm import unit

_INSTALLED = False


class _OEHasAtomName:
    def __init__(self, name: str) -> None:
        self.name = name

    def __call__(self, atom: _OEAtom) -> bool:
        return atom.GetName() == self.name


class _OEAtomPredicate:
    def __init__(self, func) -> None:
        self._func = func

    def __call__(self, atom: _OEAtom) -> bool:
        return self._func(atom)


class _OEAndAtom:
    def __init__(self, *predicates) -> None:
        self._predicates = predicates

    def __call__(self, atom: _OEAtom) -> bool:
        return all(predicate(atom) for predicate in self._predicates)


class _OENotBond:
    def __init__(self, predicate) -> None:
        self._predicate = predicate

    def __call__(self, bond: _OEBond) -> bool:
        return not self._predicate(bond)


class _OEIsRotor:
    def __call__(self, bond: _OEBond) -> bool:
        return bond.IsRotor()


class _OEIsHeavy:
    def __call__(self, atom: _OEAtom) -> bool:
        return atom.GetAtomicNum() > 1


class _OEIsAromaticAtom:
    def __call__(self, atom: _OEAtom) -> bool:
        return atom.IsAromatic()


class _OEMatchPairAtom:
    def __init__(self, pattern: _OEAtom, target: _OEAtom) -> None:
        self.pattern = pattern
        self.target = target


class _OEMatch:
    def __init__(self, pairs: list[_OEMatchPairAtom]) -> None:
        self._pairs = pairs

    def GetAtoms(self):
        return self._pairs


class _OEMCSMaxBondsCompleteCycles:
    pass


class _OEMCSSearch:
    def __init__(self, mcs_type=None) -> None:
        self._pattern: _OEMol | None = None
        self._constraints: list[tuple[str, str]] = []

    def Init(self, pattern: _OEMol, atomexpr, bondexpr) -> None:
        self._pattern = pattern

    def AddConstraint(self, match_pair: _OEMatchPairAtom) -> bool:
        self._constraints.append((match_pair.pattern.GetName(), match_pair.target.GetName()))
        return True

    def SetMCSFunc(self, func) -> None:
        return None

    def GetPattern(self) -> _OEMol:
        assert self._pattern is not None
        return self._pattern

    def Match(self, target: _OEMol, unique: bool = True):
        from rdkit import Chem
        from rdkit.Chem import rdFMCS

        assert self._pattern is not None
        rd_pattern = self._pattern._rdmol
        rd_target = target._rdmol
        result = rdFMCS.FindMCS(
            [rd_pattern, rd_target],
            atomCompare=rdFMCS.AtomCompare.Elements,
            bondCompare=rdFMCS.BondCompare.Order,
            ringMatchesRing=True,
            completeRingsOnly=False,
        )
        if result.canceled or not result.smartsString:
            return []

        query = Chem.MolFromSmarts(result.smartsString)
        pattern_matches = rd_pattern.GetSubstructMatches(query, uniquify=True)
        target_matches = rd_target.GetSubstructMatches(query, uniquify=True)
        if not pattern_matches or not target_matches:
            return []

        matches: list[_OEMatch] = []
        for pattern_match in pattern_matches:
            for target_match in target_matches:
                pairs: list[_OEMatchPairAtom] = []
                ok = True
                for pattern_idx, target_idx in zip(pattern_match, target_match):
                    pattern_atom = self._pattern.GetAtomByIdx(int(pattern_idx))
                    target_atom = target.GetAtomByIdx(int(target_idx))
                    pairs.append(_OEMatchPairAtom(pattern_atom, target_atom))
                for name_a, name_b in self._constraints:
                    if not any(
                        pair.pattern.GetName() == name_a and pair.target.GetName() == name_b
                        for pair in pairs
                    ):
                        ok = False
                        break
                if ok:
                    matches.append(_OEMatch(pairs))
        if unique and matches:
            return [matches[0]]
        return matches


class _OEBond:
    def __init__(self, mol: _OEMol, rd_bond, idx: int) -> None:
        self._mol = mol
        self._rd_bond = rd_bond
        self._idx = idx

    def GetBgn(self) -> _OEAtom:
        return self._mol.GetAtomByIdx(self._rd_bond.GetBeginAtomIdx())

    def GetEnd(self) -> _OEAtom:
        return self._mol.GetAtomByIdx(self._rd_bond.GetEndAtomIdx())

    def GetBgnIdx(self) -> int:
        return self.GetBgn().GetIdx()

    def GetEndIdx(self) -> int:
        return self.GetEnd().GetIdx()

    def GetNbr(self, atom: _OEAtom) -> _OEAtom:
        other = self._rd_bond.GetOtherAtom(atom._rd_atom)
        return self._mol.GetAtomByIdx(other.GetIdx())

    def IsRotor(self) -> bool:
        return self._rd_bond.GetBondTypeAsDouble() == 1.0


class _OETorsion:
    def __init__(self, atoms: tuple[_OEAtom, _OEAtom, _OEAtom, _OEAtom], radians: float) -> None:
        self.a, self.b, self.c, self.d = atoms
        self.radians = radians


class _OEAtom:
    def __init__(self, mol: _OEMol, rd_atom, idx: int) -> None:
        self._mol = mol
        self._rd_atom = rd_atom
        self._idx = idx
        self._data: dict[str, object] = {}
        self._name_override: str | None = None

    def GetIdx(self) -> int:
        return self._idx

    def GetName(self) -> str:
        if self._name_override is not None:
            return self._name_override
        info = self._rd_atom.GetPDBResidueInfo()
        if info is not None and info.GetName():
            return info.GetName().strip()
        if self._rd_atom.HasProp("_TriposAtomName"):
            return self._rd_atom.GetProp("_TriposAtomName")
        return self._rd_atom.GetSymbol()

    def SetName(self, name: str) -> None:
        self._name_override = name
        self._rd_atom.SetProp("_TriposAtomName", name)
        info = self._rd_atom.GetPDBResidueInfo()
        if info is not None:
            info.SetName(name.ljust(4)[:4])

    def GetAtomicNum(self) -> int:
        return self._rd_atom.GetAtomicNum()

    def IsChiral(self) -> bool:
        from rdkit import Chem

        return self._rd_atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED

    def HasStereoSpecified(self) -> bool:
        from rdkit import Chem

        return self._rd_atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED

    def IsAromatic(self) -> bool:
        return self._rd_atom.GetIsAromatic()

    def GetBonds(self):
        for bond in self._rd_atom.GetBonds():
            yield self._mol._bond_by_rd_id[bond.GetIdx()]

    def GetAtoms(self):
        for nbr in self._rd_atom.GetNeighbors():
            yield self._mol.GetAtomByIdx(nbr.GetIdx())

    def GetData(self, key: str):
        return self._data.get(key)

    def SetData(self, key: str, value) -> None:
        self._data[key] = value

    def AddData(self, key: str, value) -> None:
        self._data[key] = value


class _OEMol:
    OEMol = None  # patched after class definition

    def __init__(self, rdmol) -> None:
        from rdkit import Chem

        self._rdmol = Chem.Mol(rdmol)
        if self._rdmol.GetNumConformers() == 0 and rdmol.GetNumConformers() > 0:
            self._rdmol.AddConformer(rdmol.GetConformer(), assignId=True)
        self._atoms = [_OEAtom(self, atom, idx) for idx, atom in enumerate(self._rdmol.GetAtoms())]
        self._bond_by_rd_id = {}
        self._bonds = []
        for rd_bond in self._rdmol.GetBonds():
            bond = _OEBond(self, rd_bond, len(self._bonds))
            self._bonds.append(bond)
            self._bond_by_rd_id[rd_bond.GetIdx()] = bond
        self._title = self._rdmol.GetProp("_Name") if self._rdmol.HasProp("_Name") else "MOL"

    @classmethod
    def from_file(cls, path: str | Path, *, allow_undefined_stereo: bool = False) -> _OEMol:
        path = Path(path)
        from rdkit import Chem
        from openff.toolkit import Molecule as OFFMolecule

        if path.suffix.lower() == ".sdf":
            off = OFFMolecule.from_file(str(path), allow_undefined_stereo=allow_undefined_stereo)
            rdmol = off.to_rdkit()
            if off.n_conformers:
                conf = off.conformers[0].m_as("angstrom")
                if rdmol.GetNumConformers() == 0:
                    from rdkit.Geometry import Point3D

                    rd_conf = Chem.Conformer(rdmol.GetNumAtoms())
                    for idx, xyz in enumerate(conf):
                        rd_conf.SetAtomPosition(idx, Point3D(float(xyz[0]), float(xyz[1]), float(xyz[2])))
                    rdmol.AddConformer(rd_conf, assignId=True)
            return cls(rdmol)

        if path.suffix.lower() == ".pdb":
            rdmol = Chem.MolFromPDBFile(str(path), removeHs=False, sanitize=False)
            if rdmol is None:
                raise ValueError(f"RDKit failed to read PDB: {path}")
            Chem.SanitizeMol(rdmol, catchErrors=True)
            if rdmol.GetNumConformers() == 0:
                raise ValueError(f"PDB has no coordinates: {path}")
            return cls(rdmol)

        raise ValueError(f"Unsupported molecule file for OpenEye shim: {path}")

    def NumAtoms(self) -> int:
        return len(self._atoms)

    def GetAtoms(self, predicate=None):
        atoms = self._atoms
        if predicate is None:
            return iter(atoms)
        return (atom for atom in atoms if predicate(atom))

    def GetBonds(self):
        return iter(self._bonds)

    def GetAtom(self, predicate) -> _OEAtom:
        if isinstance(predicate, _OEHasAtomName):
            predicate = predicate
        for atom in self._atoms:
            if callable(predicate) and predicate(atom):
                return atom
        raise KeyError(f"No atom matched predicate {predicate!r}")

    def GetAtomByIdx(self, idx: int) -> _OEAtom:
        return self._atoms[idx]

    def GetTitle(self) -> str:
        return self._title

    def SetTitle(self, title: str) -> None:
        self._title = title

    def GetCoords(self) -> dict[int, tuple[float, float, float]]:
        conf = self._rdmol.GetConformer()
        return {
            idx: (
                float(conf.GetAtomPosition(idx).x),
                float(conf.GetAtomPosition(idx).y),
                float(conf.GetAtomPosition(idx).z),
            )
            for idx in range(self._rdmol.GetNumAtoms())
        }

    def SetCoords(self, coords: dict[int, tuple[float, float, float]]) -> None:
        conf = self._rdmol.GetConformer()
        for idx, xyz in coords.items():
            conf.SetAtomPosition(int(idx), xyz)


class _OEOmega:
    def SetMaxConfs(self, value: int) -> None:
        return None

    def SetStrictStereo(self, value: bool) -> None:
        return None

    def SetIncludeInput(self, value: bool) -> None:
        return None

    def __call__(self, molecule: _OEMol) -> bool:
        return True


def _build_oeomega_module() -> types.ModuleType:
    oeomega = types.ModuleType("oeomega")
    oeomega.OEOmega = _OEOmega
    oeomega.OEOmegaIsLicensed = lambda: False
    return oeomega


def _build_optional_openeye_module(name: str, license_func: str) -> types.ModuleType:
    module = types.ModuleType(name)
    setattr(module, license_func, lambda: False)
    return module


def _oemol_to_openmm_topology(molecule: _OEMol):
    from openmm.app import Element, Topology

    topology = Topology()
    chain = topology.addChain()
    residue = topology.addResidue(molecule.GetTitle() or "MOL", chain)
    for atom in molecule.GetAtoms():
        topology.addAtom(atom.GetName(), Element.getByAtomicNumber(atom.GetAtomicNum()), residue)
    atoms_by_name = {atom.name: atom for atom in topology.atoms()}
    for bond in molecule.GetBonds():
        a = bond.GetBgn().GetName()
        b = bond.GetEnd().GetName()
        topology.addBond(atoms_by_name[a], atoms_by_name[b])
    return topology


def _extract_positions(molecule: _OEMol, units=unit.angstrom):
    coords = molecule.GetCoords()
    positions = unit.Quantity(np.zeros((molecule.NumAtoms(), 3), dtype=np.float32), units)
    for idx in range(molecule.NumAtoms()):
        positions[idx, :] = unit.Quantity(coords[idx], units)
    return positions


def _describe_oemol(mol: _OEMol) -> str:
    lines = ["ATOMS:"]
    for atom in mol.GetAtoms():
        lines.append(f"{atom.GetIdx():8d} {atom.GetName():>5s} {atom.GetAtomicNum():5d}")
    lines.append("BONDS:")
    for bond in mol.GetBonds():
        lines.append(f"{bond.GetBgn().GetIdx():8d} {bond.GetEnd().GetIdx():8d}")
    return "\n".join(lines) + "\n"


def _create_oemol_from_sdf(
    sdf_filename,
    index: int = 0,
    add_hydrogens: bool = True,
    allow_undefined_stereo: bool = False,
) -> _OEMol:
    from rdkit import Chem

    path = Path(str(sdf_filename))
    if path.suffix.lower() == ".pdb":
        mol = _OEMol.from_file(path)
    else:
        supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=False)
        molecules = [entry for entry in supplier if entry is not None]
        if not molecules:
            raise ValueError(f"No molecules read from {path}")
        mol = _OEMol(molecules[index])
        if add_hydrogens:
            mol._rdmol = Chem.AddHs(mol._rdmol, addCoords=True)
            mol = _OEMol(mol._rdmol)
    names = [atom.GetName() for atom in mol.GetAtoms()]
    if len(names) != len(set(names)):
        for idx, atom in enumerate(mol.GetAtoms()):
            atom.SetName(f"{atom.GetName()}{idx}")
    return mol


def _oe_get_torsions(molecule: _OEMol, predicate) -> list[_OETorsion]:
    from rdkit.Chem import rdMolTransforms

    torsions: list[_OETorsion] = []
    rd_mol = molecule._rdmol
    conf = rd_mol.GetConformer()
    for rd_torsion in itertools.permutations(range(rd_mol.GetNumAtoms()), 4):
        a, b, c, d = rd_torsion
        if not (
            rd_mol.GetBondBetweenAtoms(a, b)
            and rd_mol.GetBondBetweenAtoms(b, c)
            and rd_mol.GetBondBetweenAtoms(c, d)
        ):
            continue
        bond_bc = molecule._bond_by_rd_id[rd_mol.GetBondBetweenAtoms(b, c).GetIdx()]
        if not predicate(bond_bc):
            continue
        if rd_mol.GetAtomWithIdx(a).GetAtomicNum() == 1 or rd_mol.GetAtomWithIdx(d).GetAtomicNum() == 1:
            continue
        radians = rdMolTransforms.GetDihedralRad(conf, a, b, c, d)
        atoms = tuple(molecule.GetAtomByIdx(i) for i in (a, b, c, d))
        torsions.append(_OETorsion(atoms, float(radians)))
    return torsions


def _oe_get_angle(molecule: _OEMol, a: _OEAtom, b: _OEAtom, c: _OEAtom) -> float:
    from rdkit.Chem import rdMolTransforms

    conf = molecule._rdmol.GetConformer()
    return float(
        rdMolTransforms.GetAngleRad(
            conf,
            a.GetIdx(),
            b.GetIdx(),
            c.GetIdx(),
        )
    )


def _oe_perceive_cip_stereo(molecule: _OEMol, atom: _OEAtom) -> int:
    return 0


def _oe_set_cip_stereo(molecule: _OEMol, atom: _OEAtom, stereo: int) -> None:
    return None


def _build_oechem_module() -> types.ModuleType:
    oechem = types.ModuleType("oechem")
    oechem.OEChemIsLicensed = lambda: False
    oechem.OEMol = _OEMol
    oechem.OEGraphMol = lambda mol: _OEMol(mol._rdmol)
    oechem.OEHasAtomName = _OEHasAtomName
    oechem.OEMCSSearch = _OEMCSSearch
    oechem.OEMCSType_Exhaustive = object()
    oechem.OEExprOpts_Aromaticity = 1
    oechem.OEExprOpts_RingMember = 2
    oechem.OEExprOpts_Degree = 4
    oechem.OEExprOpts_AtomicNumber = 8
    oechem.OEMatchPairAtom = _OEMatchPairAtom
    oechem.OEMCSMaxBondsCompleteCycles = _OEMCSMaxBondsCompleteCycles
    oechem.OEIsRotor = _OEIsRotor
    oechem.OENotBond = _OENotBond
    oechem.OEGetTorsions = _oe_get_torsions
    oechem.OEPerceiveCIPStereo = _oe_perceive_cip_stereo
    oechem.OESetCIPStereo = _oe_set_cip_stereo
    oechem.OEIsAromaticAtom = _OEIsAromaticAtom
    oechem.OEIsHeavy = _OEIsHeavy
    oechem.OEAndAtom = _OEAndAtom
    oechem.OEGetAngle = _oe_get_angle
    return oechem


def install_openeye_shim() -> None:
    """Register fake OpenEye modules and patch Perses/OpenFF entry points."""
    global _INSTALLED
    if _INSTALLED:
        return

    oechem = _build_oechem_module()
    oeomega = _build_oeomega_module()
    oequacpac = _build_optional_openeye_module("oequacpac", "OEQuacPacIsLicensed")
    oeiupac = _build_optional_openeye_module("oeiupac", "OEIUPACIsLicensed")

    openeye = types.ModuleType("openeye")
    openeye.__version__ = "shim.0.0"
    openeye.oechem = oechem
    openeye.oeomega = oeomega
    openeye.oequacpac = oequacpac
    openeye.oeiupac = oeiupac
    sys.modules["openeye"] = openeye
    sys.modules["openeye.oechem"] = oechem
    sys.modules["openeye.oeomega"] = oeomega
    sys.modules["openeye.oequacpac"] = oequacpac
    sys.modules["openeye.oeiupac"] = oeiupac

    import openmoltools.forcefield_generators as ff_generators
    import perses.utils.openeye as perses_openeye
    from openff.toolkit.topology import Molecule

    ff_generators.generateTopologyFromOEMol = lambda molecule: _oemol_to_openmm_topology(molecule)
    perses_openeye.createOEMolFromSDF = _create_oemol_from_sdf
    perses_openeye.extractPositionsFromOEMol = _extract_positions
    perses_openeye.describe_oemol = _describe_oemol

    original_from_openeye = Molecule.from_openeye.__func__

    @classmethod
    def from_openeye(cls, oemol, allow_undefined_stereo=False):
        if isinstance(oemol, _OEMol):
            return cls.from_rdkit(oemol._rdmol, allow_undefined_stereo=allow_undefined_stereo)
        return original_from_openeye(oemol, allow_undefined_stereo=allow_undefined_stereo)

    Molecule.from_openeye = from_openeye
    _INSTALLED = True


def openeye_shim_installed() -> bool:
    return _INSTALLED
