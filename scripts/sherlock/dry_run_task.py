#!/usr/bin/env python3
"""Dry-run validation for one MD manifest task (no simulation stepping).

This checks whether a task can be initialized on the target runtime by validating:
- manifest/task lookup
- prepared topology + system readability
- topology/system particle count consistency
- OpenMM simulation context creation (heating + production systems)
- stripped analysis topology + DCD reporter construction
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
import tempfile
from pathlib import Path


def _setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s %(levelname)s %(message)s")


def _safe_simulation(app, topology, system, integrator, platform, properties, allow_fallback: bool):
    try:
        return app.Simulation(topology, system, integrator, platform, properties)
    except Exception:
        if not allow_fallback:
            raise
        logging.warning("Primary platform initialization failed; falling back to default platform.")
        return app.Simulation(topology, system, integrator)


def main(argv: list[str] | None = None) -> int:
    # Ensure repository root is importable when script is invoked by path.
    repo_root = Path(__file__).resolve().parents[2]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

    parser = argparse.ArgumentParser(description="Dry-run one MD task without stepping dynamics.")
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--task-id", type=int, required=True)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args(argv)

    _setup_logging(args.verbose)

    from src.cluster.manifest import get_task_by_id
    from src.openmm.md_protocol import (
        MDProtocolConfig,
        _StrippedDCDReporter,
        _solute_atom_indices,
        _write_stripped_topology_pdb,
    )
    from src.openmm.platform import get_platform
    from src.openmm.require import require_module

    task = get_task_by_id(args.manifest, args.task_id)
    logging.info("Loaded task_id=%d mutation=%s replicate=%d", task.task_id, task.mutation, task.replicate)

    if not task.prepared_topology_pdb or not task.prepared_system_xml:
        raise ValueError("Task is missing prepared_topology_pdb/prepared_system_xml")

    topology_pdb = Path(task.prepared_topology_pdb)
    system_xml = Path(task.prepared_system_xml)
    if not topology_pdb.exists():
        raise FileNotFoundError(f"Missing topology PDB: {topology_pdb}")
    if not system_xml.exists():
        raise FileNotFoundError(f"Missing system XML: {system_xml}")

    app = require_module("openmm.app")
    openmm = require_module("openmm")
    unit = require_module("openmm.unit")

    cfg = MDProtocolConfig()

    with open(topology_pdb, "r") as handle:
        pdb = app.PDBFile(handle)
    system = openmm.XmlSerializer.deserialize(system_xml.read_text())

    n_top_atoms = sum(1 for _ in pdb.topology.atoms())
    n_sys_particles = system.getNumParticles()
    if n_top_atoms != n_sys_particles:
        raise RuntimeError(
            f"Topology/system mismatch: topology atoms={n_top_atoms}, system particles={n_sys_particles}"
        )
    logging.info("Topology/system check OK: %d atoms/particles", n_top_atoms)

    platform, properties = get_platform()
    logging.info(
        "Platform selection: available requested=%s CUDA_VISIBLE_DEVICES=%s",
        platform.getName(),
        os.environ.get("CUDA_VISIBLE_DEVICES"),
    )

    allow_fallback = str(os.environ.get("OPENMM_ALLOW_FALLBACK", "0")).strip().lower() in {
        "1",
        "true",
        "yes",
    }

    # Heating-style context check.
    heat_integrator = openmm.LangevinMiddleIntegrator(
        cfg.temperature_target_k * unit.kelvin,
        1.0 / unit.picosecond,
        cfg.timestep_fs * unit.femtoseconds,
    )
    sim = _safe_simulation(
        app, pdb.topology, system, heat_integrator, platform, properties, allow_fallback
    )
    sim.context.setPositions(pdb.positions)
    _ = sim.context.getState(getEnergy=True)
    logging.info("Heating-style context check OK")
    del sim

    # Production-style context check (with barostat).
    prod_system = openmm.XmlSerializer.deserialize(system_xml.read_text())
    prod_system.addForce(
        openmm.MonteCarloBarostat(
            cfg.pressure_bar * unit.bar,
            cfg.temperature_target_k * unit.kelvin,
            25,
        )
    )
    prod_integrator = openmm.LangevinMiddleIntegrator(
        cfg.temperature_target_k * unit.kelvin,
        1.0 / unit.picosecond,
        cfg.timestep_fs * unit.femtoseconds,
    )
    prod = _safe_simulation(
        app, pdb.topology, prod_system, prod_integrator, platform, properties, allow_fallback
    )
    prod.context.setPositions(pdb.positions)
    _ = prod.context.getState(getEnergy=True)
    logging.info("Production-style context check OK")
    del prod

    # Analysis-topology + stripped DCD reporter check.
    solute_idx = _solute_atom_indices(pdb.topology)
    logging.info("Found %d solute atoms for analysis reporter", len(solute_idx))
    with tempfile.TemporaryDirectory(prefix="nnrti_dry_run_") as tmpdir:
        tmp = Path(tmpdir)
        analysis_top = tmp / "analysis_topology.pdb"
        analysis_dcd = tmp / "analysis.dcd"
        _write_stripped_topology_pdb(pdb.topology, pdb.positions, analysis_top)
        with open(analysis_top, "r") as handle:
            stripped_topology = app.PDBFile(handle).topology
        reporter = _StrippedDCDReporter(
            file_path=analysis_dcd,
            stripped_topology=stripped_topology,
            timestep_ps=cfg.timestep_fs / 1000.0,
            interval=max(1, int(cfg.report_interval_steps)),
            atom_indices=solute_idx,
        )
        reporter.close()
    logging.info("Analysis reporter check OK")

    logging.info("Dry-run PASSED for task_id=%d (%s rep%d)", task.task_id, task.mutation, task.replicate)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
