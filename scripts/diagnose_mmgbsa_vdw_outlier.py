#!/usr/bin/env python3
"""Diagnose per-frame MM/GBSA vdW outliers for selected trajectories."""
from __future__ import annotations

import argparse
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd

from src.analysis.cli.compute_mmgbsa_safe import (
    _infer_rep_dir,
    _infer_total_steps,
    _nonempty_path,
    _resolve_local_path,
)
from src.analysis.result_collector import collect_md_results
from src.md.openmm.ligand import build_forcefield, load_ligand_molecule
from src.md.openmm.mmgbsa import (
    _apply_h_relax,
    _build_component_system,
    _build_h_relax_context,
    _energy_of,
    _extract_ligand_indices,
    _extract_receptor_indices,
    _make_context,
    _make_subtopology,
    _select_snapshot_indices,
)
from src.md.openmm.require import require_module


def _atom_label(topology_atoms: list, index: int) -> str:
    atom = topology_atoms[int(index)]
    residue = atom.residue
    return f"{int(index)}:{residue.chain.id}:{residue.name}{residue.id}:{atom.name}"


def _closest_pair(pos_nm: np.ndarray, left: np.ndarray, right: np.ndarray) -> tuple[float, int, int]:
    delta = pos_nm[left[:, None], :] - pos_nm[right[None, :], :]
    dist = np.sqrt(np.sum(delta * delta, axis=2))
    flat = int(np.nanargmin(dist))
    i, j = np.unravel_index(flat, dist.shape)
    return float(dist[i, j]), int(left[i]), int(right[j])


def _resolve_run(args: argparse.Namespace) -> tuple[pd.Series, Path, Path, Path, Path, float | None]:
    md_df = collect_md_results(args.manifest, args.results_dir)
    sub = md_df[
        md_df["mutation"].astype(str).eq(str(args.mutation))
        & (pd.to_numeric(md_df["replicate"], errors="coerce").astype("Int64") == int(args.replicate))
    ]
    if sub.empty:
        raise ValueError(f"No run found for {args.mutation} rep{args.replicate}")
    row = sub.iloc[0]
    rep = int(row["replicate"])
    safe = str(row["safe_label"])
    rep_dir = _infer_rep_dir(row)

    min_pdb = _resolve_local_path(
        _nonempty_path(row.get("minimized_pdb")),
        rep_dir / f"{safe}_minimized_rep{rep:02d}.pdb",
    )
    dcd = _resolve_local_path(
        _nonempty_path(row.get("analysis_dcd")),
        rep_dir / f"{safe}_rep{rep:02d}_analysis.dcd",
    )
    analysis_topo = _resolve_local_path(
        _nonempty_path(row.get("analysis_topology_pdb")),
        rep_dir / f"{safe}_rep{rep:02d}_analysis_topology.pdb",
    )
    ligand_sdf = _resolve_local_path(_nonempty_path(row.get("ligand_sdf")))
    if min_pdb is None or dcd is None or analysis_topo is None or ligand_sdf is None:
        raise FileNotFoundError("Missing MM/GBSA inputs")

    total_steps = _infer_total_steps(row, rep_dir, safe, rep)
    total_time_ns = float(total_steps) * 2.0 / 1_000_000.0 if total_steps else None
    return row, min_pdb, dcd, analysis_topo, ligand_sdf, total_time_ns


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=Path("manifests/md_manifest.csv"))
    parser.add_argument("--results-dir", type=Path, default=Path("results"))
    parser.add_argument("--mutation", type=str, default="Y188L")
    parser.add_argument("--replicate", type=int, default=1)
    parser.add_argument("--ligand-resname", type=str, default="2KW")
    parser.add_argument("--sample-window-ns", type=float, default=0.0)
    parser.add_argument("--sample-last-frames", type=int, default=20)
    parser.add_argument("--snapshots", type=int, default=100)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/analysis/binding_energy/vdw_diagnostics/y188l_rep1_vdw_frames.csv"),
    )
    args = parser.parse_args()

    _row, _min_pdb, dcd, topology_pdb, ligand_sdf, total_time_ns = _resolve_run(args)
    app = require_module("openmm.app")
    with open(topology_pdb) as handle:
        pdb = app.PDBFile(handle)
    topology = pdb.topology
    topology_atoms = list(topology.atoms())
    ligand = load_ligand_molecule(ligand_sdf)
    forcefield = build_forcefield([ligand])

    receptor_idx = _extract_receptor_indices(topology, args.ligand_resname)
    ligand_idx = _extract_ligand_indices(topology, args.ligand_resname)
    receptor_top = _make_subtopology(topology, pdb.positions, args.ligand_resname, invert=True)
    ligand_top = _make_subtopology(topology, pdb.positions, args.ligand_resname, invert=False)

    systems = {
        "complex_vdw": _build_component_system(topology, forcefield, "vdw"),
        "receptor_vdw": _build_component_system(receptor_top, forcefield, "vdw"),
        "ligand_vdw": _build_component_system(ligand_top, forcefield, "vdw"),
    }
    contexts = {name: _make_context(system) for name, system in systems.items()}
    h_relax_ctx, h_relax_force, h_relax_heavy = _build_h_relax_context(topology, forcefield)

    universe = mda.Universe(str(topology_pdb), str(dcd))
    window_ns = float(args.sample_window_ns) if float(args.sample_window_ns) > 0 else None
    last_frames = int(args.sample_last_frames) if int(args.sample_last_frames) > 0 else None
    frame_indices = _select_snapshot_indices(
        n_frames=len(universe.trajectory),
        discard_fraction=0.25,
        n_snapshots=int(args.snapshots),
        dt_ps=getattr(universe.trajectory, "dt", None),
        sample_window_ns=window_ns,
        total_time_ns=total_time_ns,
        sample_last_frames=last_frames,
    )

    rows = []
    for frame in frame_indices:
        universe.trajectory[int(frame)]
        pos_before = np.asarray(universe.atoms.positions, dtype=float) / 10.0
        min_before, rec_before, lig_before = _closest_pair(pos_before, receptor_idx, ligand_idx)
        pos_after = _apply_h_relax(h_relax_ctx, h_relax_force, h_relax_heavy, pos_before)
        min_after, rec_after, lig_after = _closest_pair(pos_after, receptor_idx, ligand_idx)
        rec_pos = pos_after[receptor_idx]
        lig_pos = pos_after[ligand_idx]

        complex_vdw = _energy_of(contexts["complex_vdw"], pos_after, 1)
        receptor_vdw = _energy_of(contexts["receptor_vdw"], rec_pos, 1)
        ligand_vdw = _energy_of(contexts["ligand_vdw"], lig_pos, 1)
        binding_vdw = complex_vdw - receptor_vdw - ligand_vdw

        time_ns = (
            float(frame) * float(total_time_ns) / float(len(universe.trajectory) - 1)
            if total_time_ns is not None and len(universe.trajectory) > 1
            else float("nan")
        )
        rows.append(
            {
                "mutation": args.mutation,
                "replicate": int(args.replicate),
                "frame": int(frame),
                "time_ns": time_ns,
                "complex_vdw": complex_vdw,
                "receptor_vdw": receptor_vdw,
                "ligand_vdw": ligand_vdw,
                "binding_vdw": binding_vdw,
                "closest_receptor_ligand_nm_before_h_relax": min_before,
                "closest_receptor_atom_before": _atom_label(topology_atoms, rec_before),
                "closest_ligand_atom_before": _atom_label(topology_atoms, lig_before),
                "closest_receptor_ligand_nm_after_h_relax": min_after,
                "closest_receptor_atom_after": _atom_label(topology_atoms, rec_after),
                "closest_ligand_atom_after": _atom_label(topology_atoms, lig_after),
            }
        )

    out = pd.DataFrame(rows).sort_values("frame", kind="stable")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print(out.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
