#!/usr/bin/env python
"""Doravirine alone in explicit solvent.

Reviewer 1 (comments 1 and 3) asks whether resistance acts through terms that
ddG_bind omits by construction. Two of the three legs of the binding cycle are
now covered: the bound complex (holo) and the unliganded receptor (apo). This
is the third, the free ligand in water, and it supplies two things the panel
has no reference for:

    dE_strain = <E_intra(DOR)>_bound - <E_intra(DOR)>_solution
    the free-solution distributions of tau1-tau4, against which the bound
    torsions in results/analysis/dor_torsions measure conformational restriction

DOR is identical in every genotype, so one ensemble is the common reference for
all twenty.

Everything that touches the ligand is held identical to the production runs, or
the strain term is not a difference of comparable quantities: data/ligands/dor.sdf,
Gasteiger charges, openff-2.0.0, TIP3P, 1.0 nm padding, 0.15 M NaCl, PME with a
1.0 nm cutoff, HBonds constraints, 2 fs, 25 ps of heating, NPT at 300 K / 1 bar.

Atom names are taken from the production analysis topology rather than from the
SDF, because compute_dor_torsions selects its dihedral atoms by name (C12x, C2x,
O1x, ...). The SDF and the production topology agree atom for atom in element
order, which `prep` asserts before relying on it.

Replicates are decorrelated by the thermostat: run_prepared_md calls
setVelocitiesToTemperature without a seed, so separate invocations draw
independent velocities. The 0.1 A coordinate jitter used for the production
replicates is not applied, as there is no pocket to decorrelate against.

Usage:
    python -m nnrti.cli.run_dor_solvent prep
    python -m nnrti.cli.run_dor_solvent run --replicate 1
"""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import openmm
from openmm import app, unit

from nnrti.md.openmm.ligand import build_forcefield, load_ligand_molecule
from nnrti.md.openmm.md_protocol import MDProtocolConfig, run_prepared_md
from nnrti.md.openmm.platform import get_platform

LIGAND_RESNAME = "2KW"
#: Larger than the 1.0 nm the complex is prepared with, and deliberately so.
#: DOR spans 1.21 nm, so 1.0 nm padding gives a 2.27 nm box in which the ligand
#: occupies half the edge: its periodic images sit 1.06 nm apart, barely outside
#: the 1.0 nm cutoff, and half the box exceeds the cutoff by only 12%, which an
#: NPT volume fluctuation in a 321-water box can erase. Both artifacts bias the
#: torsional populations this ensemble exists to measure. The force field,
#: charges, water model, cutoff and constraints are unchanged.
SOLVENT_PADDING_NM = 1.5
DEFAULT_LIGAND_SDF = Path("data/ligands/dor.sdf")
DEFAULT_ROOT = Path("results/md_runs/dor_solvent")
#: Production topology the atom names are copied from, so that the torsion
#: analysis finds the same names in the solution runs as in the bound runs.
NAME_REFERENCE_PDB = Path("results/md_runs/wt/rep_02/wt_rep02_analysis_topology.pdb")

logger = logging.getLogger(__name__)


def _reference_atom_names(pdb_path: Path, resname: str = LIGAND_RESNAME) -> list[str]:
    """Atom names of the ligand, in file order, from a production topology PDB."""
    names: list[str] = []
    with pdb_path.open() as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")) and line[17:20].strip() == resname:
                names.append(line[12:16].strip())
    if not names:
        raise ValueError(f"No {resname} atoms found in {pdb_path}")
    return names


def _ligand_topology(molecule, atom_names: list[str]):
    """OpenMM topology for the ligand alone, with production atom names.

    Built from the OpenFF molecule rather than parsed from a PDB so the bond
    graph is explicit: the SMIRNOFF template generator matches on connectivity,
    and a HETATM-only PDB without CONECT records would not carry it.
    """
    off_top = molecule.to_topology()
    topology = off_top.to_openmm()

    elements = [a.element.symbol for a in topology.atoms()]
    sdf_elements = [a.symbol for a in molecule.atoms]
    if elements != sdf_elements:
        raise AssertionError(
            "OpenMM topology atom order diverged from the SDF; refusing to assign names by position"
        )
    if len(atom_names) != len(elements):
        raise AssertionError(
            f"Name reference has {len(atom_names)} atoms, ligand has {len(elements)}"
        )

    for residue in topology.residues():
        residue.name = LIGAND_RESNAME
    for atom, name in zip(topology.atoms(), atom_names):
        atom.name = name

    positions = off_top.get_positions().to_openmm()
    return topology, positions


def prep(args: argparse.Namespace) -> int:
    cfg = MDProtocolConfig()
    out_dir = args.root / "prepared"
    out_dir.mkdir(parents=True, exist_ok=True)
    topology_pdb = out_dir / "dor_solvent_topology.pdb"
    system_xml = out_dir / "dor_solvent_system.xml"

    names = _reference_atom_names(args.name_reference)
    logger.info("Atom names taken from %s (%d atoms)", args.name_reference, len(names))

    molecule = load_ligand_molecule(args.ligand_sdf)
    forcefield = build_forcefield([molecule])
    lig_top, lig_pos = _ligand_topology(molecule, names)

    modeller = app.Modeller(lig_top, lig_pos)
    modeller.addSolvent(
        forcefield,
        model="tip3p",
        padding=args.padding_nm * unit.nanometer,
        ionicStrength=cfg.ionic_strength_molar * unit.molar,
    )
    logger.info(
        "Solvated system: %d atoms (padding %.2f nm)",
        modeller.topology.getNumAtoms(),
        args.padding_nm,
    )

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=app.HBonds,
    )

    platform, properties = get_platform()
    integrator = openmm.LangevinMiddleIntegrator(
        cfg.temperature_target_k * unit.kelvin,
        1.0 / unit.picosecond,
        cfg.timestep_fs * unit.femtoseconds,
    )
    simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy(maxIterations=500)
    positions = simulation.context.getState(getPositions=True).getPositions()
    del simulation  # Release the GPU context before any later simulation.

    with topology_pdb.open("w") as handle:
        app.PDBFile.writeFile(modeller.topology, positions, handle)
    system_xml.write_text(openmm.XmlSerializer.serialize(system))

    logger.info("Wrote %s and %s", topology_pdb, system_xml)
    return 0


def run(args: argparse.Namespace) -> int:
    prepared_dir = args.root / "prepared"
    topology_pdb = prepared_dir / "dor_solvent_topology.pdb"
    system_xml = prepared_dir / "dor_solvent_system.xml"
    for path in (topology_pdb, system_xml):
        if not path.exists():
            raise FileNotFoundError(f"{path} missing; run `prep` first")

    label = f"dor_solvent_rep{args.replicate:02d}"
    out_dir = args.root / f"rep_{args.replicate:02d}"
    out_dir.mkdir(parents=True, exist_ok=True)

    timestep_fs = 2.0
    production_steps = max(1, int(round((args.production_ns * 1_000_000.0) / timestep_fs)))
    # The system is ~50x smaller than the solvated complex, so frames are cheap:
    # take ~2000 of them for torsion statistics rather than the 200 the panel uses.
    analysis_interval = max(1, production_steps // args.analysis_frames)

    cfg = MDProtocolConfig(
        production_ns=args.production_ns,
        report_interval_steps=args.report_interval,
        analysis_report_interval_steps=analysis_interval,
    )

    result = run_prepared_md(
        prepared_topology_pdb=topology_pdb,
        prepared_system_xml=system_xml,
        final_pdb_path=out_dir / f"{label}_md_final.pdb",
        state_csv_path=out_dir / f"{label}_md_state.csv",
        config=cfg,
        analysis_dcd_path=out_dir / f"{label}_analysis.dcd",
        analysis_topology_pdb_path=out_dir / f"{label}_analysis_topology.pdb",
        checkpoint_path=out_dir / f"{label}_md.chk",
        checkpoint_interval_steps=args.checkpoint_interval,
        resume_from_checkpoint=not args.no_resume,
    )

    payload = {
        "label": label,
        "replicate": args.replicate,
        "ligand_resname": LIGAND_RESNAME,
        "ligand_sdf": str(args.ligand_sdf),
        "prepared_topology_pdb": str(topology_pdb),
        "prepared_system_xml": str(system_xml),
        "md_heating_steps": result.heating_steps,
        "md_production_steps": result.production_steps,
        "md_production_steps_completed": result.production_steps_completed,
        "md_total_steps": result.total_steps,
        "analysis_interval_steps": analysis_interval,
        "resumed_from_checkpoint": result.resumed_from_checkpoint,
        "elapsed_seconds": result.elapsed_seconds,
        "status": "ok",
    }
    (out_dir / f"{label}.json").write_text(json.dumps(payload, indent=2))
    logger.info(
        "rep%02d finished %d/%d production steps in %.1f s",
        args.replicate,
        result.production_steps_completed,
        result.production_steps,
        result.elapsed_seconds,
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    p.add_argument("--ligand-sdf", type=Path, default=DEFAULT_LIGAND_SDF)
    sub = p.add_subparsers(dest="command", required=True)

    sp = sub.add_parser("prep", help="solvate and minimize DOR in water")
    sp.add_argument("--name-reference", type=Path, default=NAME_REFERENCE_PDB)
    sp.add_argument("--padding-nm", type=float, default=SOLVENT_PADDING_NM)
    sp.set_defaults(func=prep)

    sr = sub.add_parser("run", help="heating plus NPT production for one replicate")
    sr.add_argument("--replicate", type=int, required=True)
    sr.add_argument("--production-ns", type=float, default=100.0)
    sr.add_argument("--analysis-frames", type=int, default=2000)
    sr.add_argument("--report-interval", type=int, default=2000)
    sr.add_argument("--checkpoint-interval", type=int, default=5000)
    sr.add_argument("--no-resume", action="store_true")
    sr.set_defaults(func=run)

    args = p.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
