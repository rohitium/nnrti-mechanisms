#!/usr/bin/env python3
"""Compute MM/GBSA-style energy components for static PDB/mmCIF complexes."""
from __future__ import annotations

import argparse
import csv
import shlex
from pathlib import Path

import numpy as np

from src.md.openmm.ligand import build_forcefield, load_ligand_molecule
from src.md.openmm.mmgbsa import (
    _energy_of,
    _extract_ligand_indices,
    _extract_receptor_indices,
    _make_context,
    _make_subtopology,
    _radius_from_symbol,
    _subset_positions,
)
from src.md.openmm.require import require_module


def _is_hydrogen(atom) -> bool:
    if getattr(atom, "element", None) is not None:
        return atom.element.symbol.upper() == "H"
    return atom.name.upper().startswith("H")


def _chem_comp_atom_names(cif_path: Path, comp_id: str) -> list[str]:
    """Return atom names for comp_id in mmCIF _chem_comp_atom order."""
    lines = cif_path.read_text(errors="replace").splitlines()
    names: list[str] = []
    in_loop = False
    fields: list[str] = []
    for line in lines:
        stripped = line.strip()
        if stripped == "loop_":
            in_loop = True
            fields = []
            continue
        if not in_loop:
            continue
        if stripped.startswith("_chem_comp_atom."):
            fields.append(stripped)
            continue
        if fields and stripped.startswith("#"):
            if names:
                return names
            in_loop = False
            fields = []
            continue
        if fields and stripped and not stripped.startswith("_"):
            try:
                parts = shlex.split(stripped)
            except ValueError:
                continue
            if len(parts) < len(fields):
                continue
            row = dict(zip(fields, parts))
            if row.get("_chem_comp_atom.comp_id") == comp_id:
                atom_name = row.get("_chem_comp_atom.atom_id")
                if atom_name:
                    names.append(atom_name)
    return names


def _prepare_static_pdb(
    cif_path: Path,
    ligand_resname: str,
    ligand_sdf: Path,
    output_pdb: Path,
    mutation: str | None = None,
    mutation_chain: str = "A",
) -> tuple[object, object, Path]:
    """Write a protonated static complex PDB with ligand pose retained."""
    app = require_module("openmm.app")
    pdbfixer = require_module("pdbfixer")
    omm_unit = require_module("openmm.unit")
    off_unit = require_module("openff.units").unit

    ligand_mol = load_ligand_molecule(ligand_sdf)
    forcefield = build_forcefield([ligand_mol])

    with cif_path.open("r") as handle:
        fixer = pdbfixer.PDBFixer(pdbxfile=handle)
    if mutation:
        fixer.applyMutations([mutation], mutation_chain)
    for residue in fixer.topology.residues():
        if residue.name == "OMC":
            residue.name = "DC"
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)

    modeller = app.Modeller(fixer.topology, fixer.positions)

    omc_atoms = [
        atom
        for atom in modeller.topology.atoms()
        if atom.residue.name == "DC" and atom.name in {"O2'", "CM2"}
    ]
    if omc_atoms:
        modeller.delete(omc_atoms)

    # 5J2M contains MRG: N2-(3-mercaptopropyl)-2'-deoxyguanosine, covalently
    # tethered to Cys258.  The local Amber/OpenFF path has no template for this
    # crosslinked nucleotide.  For this rough static estimate, remove the single
    # crosslinked nucleotide from the receptor.
    mrg_residues = [residue for residue in modeller.topology.residues() if residue.name == "MRG"]
    if mrg_residues:
        modeller.delete(mrg_residues)

    terminal_phosphate_h = []
    for residue in modeller.topology.residues():
        atom_names = {atom.name for atom in residue.atoms()}
        if "P" in atom_names:
            terminal_phosphate_h.extend(
                atom for atom in residue.atoms() if atom.name in {"HO5'", "H5T"}
            )
    if terminal_phosphate_h:
        modeller.delete(terminal_phosphate_h)

    ligand_atoms = [
        atom for atom in modeller.topology.atoms() if atom.residue.name == ligand_resname
    ]
    if not ligand_atoms:
        raise ValueError(f"Ligand residue {ligand_resname!r} not found in {cif_path}.")

    original_heavy_by_name = {
        atom.name: modeller.positions[atom.index].value_in_unit(omm_unit.nanometer)
        for atom in ligand_atoms
        if not _is_hydrogen(atom)
    }

    ligand_residues = [
        res for res in modeller.topology.residues() if res.name == ligand_resname
    ]
    modeller.delete(ligand_residues)

    ligand_topology = ligand_mol.to_topology().to_openmm()
    for residue in ligand_topology.residues():
        residue.name = ligand_resname

    ligand_template_xyz = np.asarray(
        ligand_mol.conformers[0].to(off_unit.nanometer).magnitude,
        dtype=float,
    )
    template_heavy_mask = np.array([a.atomic_number != 1 for a in ligand_mol.atoms], dtype=bool)
    template_heavy_xyz = ligand_template_xyz[template_heavy_mask]
    ccd_heavy_names = [
        name
        for name in _chem_comp_atom_names(cif_path, ligand_resname)
        if name in original_heavy_by_name
    ]
    if len(ccd_heavy_names) == template_heavy_xyz.shape[0]:
        original_heavy_xyz = np.array(
            [original_heavy_by_name[name] for name in ccd_heavy_names],
            dtype=float,
        )
    else:
        original_heavy_xyz = np.array(list(original_heavy_by_name.values()), dtype=float)
    if template_heavy_xyz.shape[0] != original_heavy_xyz.shape[0]:
        raise ValueError(
            "Ligand heavy-atom count mismatch between CIF and SDF "
            f"({original_heavy_xyz.shape[0]} vs {template_heavy_xyz.shape[0]})."
        )

    mobile_centroid = template_heavy_xyz.mean(axis=0)
    target_centroid = original_heavy_xyz.mean(axis=0)
    m = template_heavy_xyz - mobile_centroid
    t = original_heavy_xyz - target_centroid
    u, _, vt = np.linalg.svd(m.T @ t)
    r = vt.T @ u.T
    if np.linalg.det(r) < 0:
        vt[-1, :] *= -1
        r = vt.T @ u.T
    aligned_xyz = ((ligand_template_xyz - mobile_centroid) @ r) + target_centroid
    modeller.add(ligand_topology, aligned_xyz * omm_unit.nanometer)

    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    with output_pdb.open("w") as handle:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, handle)
    return modeller.topology, modeller.positions, output_pdb


def _static_residue_templates(topology) -> dict:
    """Resolve static-prep residue template ambiguities."""
    templates = {}
    for residue in topology.residues():
        if residue.name != "CYS":
            continue
        atom_names = {atom.name for atom in residue.atoms()}
        if "HG" not in atom_names:
            templates[residue] = "CYM"
    return templates


def _build_component_system_static(topology, forcefield, mode: str):
    app = require_module("openmm.app")
    openmm = require_module("openmm")

    system = forcefield.createSystem(
        topology,
        nonbondedMethod=app.NoCutoff,
        constraints=None,
        ignoreExternalBonds=True,
        residueTemplates=_static_residue_templates(topology),
    )

    for force in system.getForces():
        force.setForceGroup(0)

    nb = None
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            nb = force
            break
    if nb is None:
        raise RuntimeError("NonbondedForce not found.")

    if mode == "elec":
        for i in range(nb.getNumParticles()):
            q, sigma, _eps = nb.getParticleParameters(i)
            nb.setParticleParameters(i, q, sigma, 0.0)
        for i in range(nb.getNumExceptions()):
            p1, p2, qprod, sigma, _eps = nb.getExceptionParameters(i)
            nb.setExceptionParameters(i, p1, p2, qprod, sigma, 0.0)
    elif mode == "vdw":
        for i in range(nb.getNumParticles()):
            _q, sigma, eps = nb.getParticleParameters(i)
            nb.setParticleParameters(i, 0.0, sigma, eps)
        for i in range(nb.getNumExceptions()):
            p1, p2, _qprod, sigma, eps = nb.getExceptionParameters(i)
            nb.setExceptionParameters(i, p1, p2, 0.0, sigma, eps)
    elif mode != "full":
        raise ValueError(f"Unknown mode: {mode}")

    nb.setForceGroup(1)
    return system


def _build_gb_system_static(topology, forcefield, include_sa: bool):
    app = require_module("openmm.app")
    openmm = require_module("openmm")

    system = forcefield.createSystem(
        topology,
        nonbondedMethod=app.NoCutoff,
        constraints=None,
        ignoreExternalBonds=True,
        residueTemplates=_static_residue_templates(topology),
    )
    for force in system.getForces():
        force.setForceGroup(0)

    nb = None
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            nb = force
            break
    if nb is None:
        raise RuntimeError("NonbondedForce not found for GB system.")

    gb = openmm.GBSAOBCForce()
    gb.setSolventDielectric(80.0)
    gb.setSoluteDielectric(2.0)
    gb.setSurfaceAreaEnergy(2.25936 if include_sa else 0.0)

    for atom in topology.atoms():
        q, _sigma, _eps = nb.getParticleParameters(atom.index)
        symbol = atom.element.symbol if atom.element is not None else atom.name[0]
        gb.addParticle(q, _radius_from_symbol(symbol), 1.0)

    gb.setForceGroup(2)
    system.addForce(gb)
    return system


def _build_full_system_static(topology, forcefield, use_cutoff: bool = False):
    app = require_module("openmm.app")
    unit = require_module("openmm.unit")

    kwargs = {}
    if use_cutoff:
        kwargs.update(
            {
                "nonbondedMethod": app.CutoffNonPeriodic,
                "nonbondedCutoff": 1.0 * unit.nanometer,
            }
        )
    else:
        kwargs["nonbondedMethod"] = app.NoCutoff

    return forcefield.createSystem(
        topology,
        constraints=None,
        ignoreExternalBonds=True,
        residueTemplates=_static_residue_templates(topology),
        **kwargs,
    )


def _minimize_static_complex(
    topology,
    positions,
    ligand_sdf: Path,
    restraint_k_kj_mol_nm2: float,
    max_iterations: int,
):
    """Minimize the prepared static complex with optional heavy-atom restraints."""
    app = require_module("openmm.app")
    openmm = require_module("openmm")
    unit = require_module("openmm.unit")

    ligand = load_ligand_molecule(ligand_sdf)
    forcefield = build_forcefield([ligand])
    system = _build_full_system_static(topology, forcefield, use_cutoff=True)

    if restraint_k_kj_mol_nm2 > 0.0:
        restraint = openmm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        restraint.addGlobalParameter("k", float(restraint_k_kj_mol_nm2))
        restraint.addPerParticleParameter("x0")
        restraint.addPerParticleParameter("y0")
        restraint.addPerParticleParameter("z0")
        for atom in topology.atoms():
            if _is_hydrogen(atom):
                continue
            pos = positions[atom.index].value_in_unit(unit.nanometer)
            restraint.addParticle(atom.index, [float(pos[0]), float(pos[1]), float(pos[2])])
        system.addForce(restraint)

    integrator = openmm.LangevinIntegrator(
        300.0 * unit.kelvin,
        1.0 / unit.picosecond,
        0.001 * unit.picoseconds,
    )
    platform_name = __import__("os").environ.get("OPENMM_PLATFORM", "CPU").strip() or "CPU"
    try:
        platform = openmm.Platform.getPlatformByName(platform_name)
        properties = {}
        if platform_name == "CPU":
            properties["Threads"] = __import__("os").environ.get("OPENMM_CPU_THREADS", "1")
        simulation = app.Simulation(topology, system, integrator, platform, properties)
    except Exception:
        simulation = app.Simulation(topology, system, integrator)

    simulation.context.setPositions(positions)
    state0 = simulation.context.getState(getEnergy=True)
    e0 = float(state0.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    openmm.LocalEnergyMinimizer.minimize(
        simulation.context,
        tolerance=10.0,
        maxIterations=max(1, int(max_iterations)),
    )
    state1 = simulation.context.getState(getPositions=True, getEnergy=True)
    e1 = float(state1.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    return state1.getPositions(), e0, e1


def _explicit_solvent_minimize(
    topology,
    positions,
    ligand_sdf: Path,
    solute_output_pdb: Path,
    solvated_output_pdb: Path,
    restraint_k_kj_mol_nm2: float,
    max_iterations: int,
):
    """Solvate with TIP3P/ions, PME-minimize, and return minimized solute positions."""
    app = require_module("openmm.app")
    openmm = require_module("openmm")
    unit = require_module("openmm.unit")

    ligand = load_ligand_molecule(ligand_sdf)
    forcefield = build_forcefield([ligand])
    modeller = app.Modeller(topology, positions)
    solute_n_atoms = len(list(topology.atoms()))

    original_create_system = forcefield.createSystem

    def _create_system_ignore_external_bonds(
        top,
        nonbondedMethod=app.NoCutoff,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=None,
        rigidWater=None,
        removeCMMotion=True,
        hydrogenMass=None,
        residueTemplates={},
        ignoreExternalBonds=False,
        switchDistance=None,
        flexibleConstraints=False,
        drudeMass=0.4 * unit.dalton,
        **args,
    ):
        return original_create_system(
            top,
            nonbondedMethod=nonbondedMethod,
            nonbondedCutoff=nonbondedCutoff,
            constraints=constraints,
            rigidWater=rigidWater,
            removeCMMotion=removeCMMotion,
            hydrogenMass=hydrogenMass,
            residueTemplates=residueTemplates or _static_residue_templates(top),
            ignoreExternalBonds=True,
            switchDistance=switchDistance,
            flexibleConstraints=flexibleConstraints,
            drudeMass=drudeMass,
            **args,
        )

    forcefield.createSystem = _create_system_ignore_external_bonds
    modeller.addSolvent(
        forcefield,
        model="tip3p",
        padding=1.0 * unit.nanometer,
        ionicStrength=0.15 * unit.molar,
        residueTemplates=_static_residue_templates(modeller.topology),
    )
    forcefield.createSystem = original_create_system

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=app.HBonds,
        ignoreExternalBonds=True,
        residueTemplates=_static_residue_templates(modeller.topology),
    )

    if restraint_k_kj_mol_nm2 > 0.0:
        restraint = openmm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        restraint.addGlobalParameter("k", float(restraint_k_kj_mol_nm2))
        restraint.addPerParticleParameter("x0")
        restraint.addPerParticleParameter("y0")
        restraint.addPerParticleParameter("z0")
        for atom in modeller.topology.atoms():
            if atom.index >= solute_n_atoms or _is_hydrogen(atom):
                continue
            pos = modeller.positions[atom.index].value_in_unit(unit.nanometer)
            restraint.addParticle(atom.index, [float(pos[0]), float(pos[1]), float(pos[2])])
        system.addForce(restraint)

    integrator = openmm.LangevinMiddleIntegrator(
        300.0 * unit.kelvin,
        1.0 / unit.picosecond,
        0.002 * unit.picoseconds,
    )
    platform_name = __import__("os").environ.get("OPENMM_PLATFORM", "CPU").strip() or "CPU"
    try:
        platform = openmm.Platform.getPlatformByName(platform_name)
        properties = {}
        if platform_name == "CPU":
            properties["Threads"] = __import__("os").environ.get("OPENMM_CPU_THREADS", "1")
        simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)
    except Exception:
        simulation = app.Simulation(modeller.topology, system, integrator)

    simulation.context.setPositions(modeller.positions)
    state0 = simulation.context.getState(getEnergy=True)
    e0 = float(state0.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    openmm.LocalEnergyMinimizer.minimize(
        simulation.context,
        tolerance=10.0,
        maxIterations=max(1, int(max_iterations)),
    )
    state1 = simulation.context.getState(getPositions=True, getEnergy=True)
    e1 = float(state1.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    solvated_positions = state1.getPositions()

    solvated_output_pdb.parent.mkdir(parents=True, exist_ok=True)
    with solvated_output_pdb.open("w") as handle:
        app.PDBFile.writeFile(modeller.topology, solvated_positions, handle)

    solute_positions = solvated_positions[:solute_n_atoms]
    with solute_output_pdb.open("w") as handle:
        app.PDBFile.writeFile(topology, solute_positions, handle)

    return solute_positions, e0, e1, solute_n_atoms, len(list(modeller.topology.atoms()))


def _compute_static_components(
    topology,
    positions,
    ligand_resname: str,
    ligand_sdf: Path,
) -> dict[str, float]:
    unit = require_module("openmm.unit")
    ligand = load_ligand_molecule(ligand_sdf)
    forcefield = build_forcefield([ligand])

    complex_top = topology
    receptor_idx = _extract_receptor_indices(complex_top, ligand_resname)
    ligand_idx = _extract_ligand_indices(complex_top, ligand_resname)
    if receptor_idx.size == 0 or ligand_idx.size == 0:
        raise ValueError("Could not identify receptor/ligand atoms.")

    receptor_top = _make_subtopology(complex_top, positions, ligand_resname, invert=True)
    ligand_top = _make_subtopology(complex_top, positions, ligand_resname, invert=False)

    systems = {
        "complex_vdw": _build_component_system_static(complex_top, forcefield, "vdw"),
        "complex_elec": _build_component_system_static(complex_top, forcefield, "elec"),
        "complex_gb_total": _build_gb_system_static(complex_top, forcefield, include_sa=True),
        "complex_gb_polar": _build_gb_system_static(complex_top, forcefield, include_sa=False),
        "receptor_vdw": _build_component_system_static(receptor_top, forcefield, "vdw"),
        "receptor_elec": _build_component_system_static(receptor_top, forcefield, "elec"),
        "receptor_gb_total": _build_gb_system_static(receptor_top, forcefield, include_sa=True),
        "receptor_gb_polar": _build_gb_system_static(receptor_top, forcefield, include_sa=False),
        "ligand_vdw": _build_component_system_static(ligand_top, forcefield, "vdw"),
        "ligand_elec": _build_component_system_static(ligand_top, forcefield, "elec"),
        "ligand_gb_total": _build_gb_system_static(ligand_top, forcefield, include_sa=True),
        "ligand_gb_polar": _build_gb_system_static(ligand_top, forcefield, include_sa=False),
    }
    contexts = {key: _make_context(system) for key, system in systems.items()}

    positions_nm = np.asarray(positions.value_in_unit(unit.nanometer), dtype=float)
    rec_nm = _subset_positions(positions_nm, receptor_idx)
    lig_nm = _subset_positions(positions_nm, ligand_idx)

    e_vdw_c = _energy_of(contexts["complex_vdw"], positions_nm, 1)
    e_vdw_r = _energy_of(contexts["receptor_vdw"], rec_nm, 1)
    e_vdw_l = _energy_of(contexts["ligand_vdw"], lig_nm, 1)
    e_elec_c = _energy_of(contexts["complex_elec"], positions_nm, 1)
    e_elec_r = _energy_of(contexts["receptor_elec"], rec_nm, 1)
    e_elec_l = _energy_of(contexts["ligand_elec"], lig_nm, 1)
    gb_t_c = _energy_of(contexts["complex_gb_total"], positions_nm, 2)
    gb_t_r = _energy_of(contexts["receptor_gb_total"], rec_nm, 2)
    gb_t_l = _energy_of(contexts["ligand_gb_total"], lig_nm, 2)
    gb_p_c = _energy_of(contexts["complex_gb_polar"], positions_nm, 2)
    gb_p_r = _energy_of(contexts["receptor_gb_polar"], rec_nm, 2)
    gb_p_l = _energy_of(contexts["ligand_gb_polar"], lig_nm, 2)

    delta_vdw = e_vdw_c - e_vdw_r - e_vdw_l
    delta_elec = e_elec_c - e_elec_r - e_elec_l
    delta_gb = gb_p_c - gb_p_r - gb_p_l
    delta_sa = (gb_t_c - gb_p_c) - (gb_t_r - gb_p_r) - (gb_t_l - gb_p_l)
    delta_total = delta_vdw + delta_elec + delta_gb + delta_sa

    return {
        "delta_g_kj_mol": float(delta_total),
        "delta_vdw_kj_mol": float(delta_vdw),
        "delta_electrostatic_kj_mol": float(delta_elec),
        "delta_gb_kj_mol": float(delta_gb),
        "delta_sa_kj_mol": float(delta_sa),
        "delta_g_kcal_mol": float(delta_total / 4.184),
        "delta_vdw_kcal_mol": float(delta_vdw / 4.184),
        "delta_electrostatic_kcal_mol": float(delta_elec / 4.184),
        "delta_gb_kcal_mol": float(delta_gb / 4.184),
        "delta_sa_kcal_mol": float(delta_sa / 4.184),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-csv", type=Path, default=Path("results/static_mmgbsa.csv"))
    parser.add_argument(
        "--prepared-dir",
        type=Path,
        default=Path("results/static_mmgbsa_prepared"),
    )
    parser.add_argument("--minimize", action="store_true", help="Minimize before MM/GBSA.")
    parser.add_argument(
        "--explicit-solvent-minimize",
        action="store_true",
        help="Solvate with TIP3P/ions and PME-minimize before MM/GBSA.",
    )
    parser.add_argument(
        "--restraint-k",
        type=float,
        default=500.0,
        help="Heavy-atom restraint k in kJ/mol/nm^2 for minimization. Use 0 for unrestrained.",
    )
    parser.add_argument("--max-iterations", type=int, default=2000)
    args = parser.parse_args()

    base_jobs = [
        {
            "structure": "7Z2G",
            "description": "RT:DNA:DOR",
            "cif": Path("data/structures/7Z2G.cif"),
            "ligand_resname": "2KW",
            "ligand_sdf": Path("data/ligands/dor.sdf"),
        },
        {
            "structure": "5J2M",
            "description": "RT:DNA:ISL-TP",
            "cif": Path("data/structures/5J2M.cif"),
            "ligand_resname": "6FN",
            "ligand_sdf": Path("data/ligands/6fn_bound.sdf"),
        },
    ]
    jobs = []
    for job in base_jobs:
        for variant, mutation in [("WT", None), ("F227C", "PHE-227-CYS")]:
            row = dict(job)
            row["variant"] = variant
            row["mutation"] = mutation
            jobs.append(row)

    rows: list[dict[str, object]] = []
    for job in jobs:
        prepared_pdb = (
            args.prepared_dir
            / f"{job['structure']}_{job['variant']}_{job['ligand_resname']}_static_prepared.pdb"
        )
        topology, positions, prepared_pdb = _prepare_static_pdb(
            cif_path=job["cif"],
            ligand_resname=str(job["ligand_resname"]),
            ligand_sdf=job["ligand_sdf"],
            output_pdb=prepared_pdb,
            mutation=job["mutation"],
        )
        min_energy_initial = float("nan")
        min_energy_final = float("nan")
        explicit_solvated_pdb = ""
        n_solute_atoms = float("nan")
        n_solvated_atoms = float("nan")
        if args.explicit_solvent_minimize:
            solvated_pdb = prepared_pdb.with_name(prepared_pdb.stem + "_solvated_minimized.pdb")
            positions, min_energy_initial, min_energy_final, n_solute_atoms, n_solvated_atoms = (
                _explicit_solvent_minimize(
                    topology=topology,
                    positions=positions,
                    ligand_sdf=job["ligand_sdf"],
                    solute_output_pdb=prepared_pdb,
                    solvated_output_pdb=solvated_pdb,
                    restraint_k_kj_mol_nm2=float(args.restraint_k),
                    max_iterations=int(args.max_iterations),
                )
            )
            explicit_solvated_pdb = str(solvated_pdb)
        elif args.minimize:
            positions, min_energy_initial, min_energy_final = _minimize_static_complex(
                topology=topology,
                positions=positions,
                ligand_sdf=job["ligand_sdf"],
                restraint_k_kj_mol_nm2=float(args.restraint_k),
                max_iterations=int(args.max_iterations),
            )
            with prepared_pdb.open("w") as handle:
                require_module("openmm.app").PDBFile.writeFile(topology, positions, handle)
        row = {
            "structure": job["structure"],
            "variant": job["variant"],
            "mutation": "F227C" if job["variant"] == "F227C" else "",
            "mutation_chain": "A" if job["variant"] == "F227C" else "",
            "minimized": bool(args.minimize or args.explicit_solvent_minimize),
            "explicit_solvent_minimized": bool(args.explicit_solvent_minimize),
            "minimization_restraint_k_kj_mol_nm2": (
                float(args.restraint_k) if (args.minimize or args.explicit_solvent_minimize) else float("nan")
            ),
            "minimization_energy_initial_kj_mol": min_energy_initial,
            "minimization_energy_final_kj_mol": min_energy_final,
            "explicit_solvent_solute_atoms": n_solute_atoms,
            "explicit_solvent_total_atoms": n_solvated_atoms,
            "description": job["description"],
            "ligand_resname": job["ligand_resname"],
            "ligand_sdf": str(job["ligand_sdf"]),
            "prepared_pdb": str(prepared_pdb),
            "solvated_minimized_pdb": explicit_solvated_pdb,
        }
        row.update(
            _compute_static_components(
                topology=topology,
                positions=positions,
                ligand_resname=str(job["ligand_resname"]),
                ligand_sdf=job["ligand_sdf"],
            )
        )
        rows.append(row)
        args.out_csv.parent.mkdir(parents=True, exist_ok=True)
        with args.out_csv.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {args.out_csv}")
    for row in rows:
        print(
            f"{row['structure']} {row['variant']} {row['ligand_resname']}: "
            f"ΔG={row['delta_g_kj_mol']:.3f} kJ/mol "
            f"({row['delta_g_kcal_mol']:.3f} kcal/mol)"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
