#!/usr/bin/env python3
"""Diagnose GB polar-solvation terms for selected MM/GBSA replicates."""
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from src.md.openmm.ligand import build_forcefield, load_ligand_molecule
from src.md.openmm.mmgbsa import (
    _apply_h_relax,
    _build_gb_system,
    _build_h_relax_context,
    _energy_of,
    _extract_ligand_indices,
    _extract_receptor_indices,
    _make_context,
    _make_subtopology,
    _select_snapshot_indices,
)
from src.md.openmm.require import require_module


@dataclass(frozen=True)
class RunInput:
    mutation: str
    safe_label: str
    replicate: int
    ligand_resname: str
    ligand_sdf: Path
    topology_pdb: Path
    trajectory_dcd: Path


def _run_input(root: Path, mutation: str, replicate: int) -> RunInput:
    safe = "wt" if mutation == "WT" else mutation.replace("+", "_")
    run_dir = root / "results" / "md_runs" / safe / f"rep_{replicate:02d}"
    json_path = run_dir / f"{safe}_rep{replicate:02d}.json"
    data = json.loads(json_path.read_text())
    return RunInput(
        mutation=mutation,
        safe_label=safe,
        replicate=replicate,
        ligand_resname=str(data["ligand_resname"]),
        ligand_sdf=root / str(data["ligand_sdf"]),
        topology_pdb=root / str(data["analysis_topology_pdb"]),
        trajectory_dcd=root / str(data["analysis_dcd"]),
    )


def _residue_atom_indices(
    topology,
    residue_number: int | None = None,
    residue_name: str | None = None,
    chain_id: str | None = None,
) -> list[int]:
    indices: list[int] = []
    for atom in topology.atoms():
        residue = atom.residue
        matches_number = residue_number is None or int(residue.id) == int(residue_number)
        matches_name = residue_name is None or residue.name.upper() == residue_name.upper()
        matches_chain = chain_id is None or str(residue.chain.id).strip() == str(chain_id).strip()
        if matches_number and matches_name and matches_chain:
            indices.append(atom.index)
    return indices


def _atom_indices_by_name(topology, residue_number: int, names: set[str], chain_id: str = "A") -> list[int]:
    indices: list[int] = []
    for atom in topology.atoms():
        if (
            int(atom.residue.id) == int(residue_number)
            and str(atom.residue.chain.id).strip() == str(chain_id).strip()
            and atom.name in names
        ):
            indices.append(atom.index)
    return indices


def _min_distance_nm(pos_nm: np.ndarray, left: list[int] | np.ndarray, right: list[int] | np.ndarray) -> float:
    left_arr = np.asarray(left, dtype=int)
    right_arr = np.asarray(right, dtype=int)
    if left_arr.size == 0 or right_arr.size == 0:
        return float("nan")
    delta = pos_nm[left_arr, None, :] - pos_nm[right_arr][None, :, :]
    distances = np.sqrt(np.sum(delta * delta, axis=2))
    return float(np.nanmin(distances))


def _evaluate_run(run: RunInput, *, n_snapshots: int, discard_fraction: float, sample_window_ns: float | None) -> pd.DataFrame:
    app = require_module("openmm.app")
    mda = require_module("MDAnalysis")

    with open(run.topology_pdb) as handle:
        solute_pdb = app.PDBFile(handle)
    topology = solute_pdb.topology
    ligand = load_ligand_molecule(run.ligand_sdf)
    forcefield = build_forcefield([ligand])

    receptor_idx = _extract_receptor_indices(topology, run.ligand_resname)
    ligand_idx = _extract_ligand_indices(topology, run.ligand_resname)
    receptor_top = _make_subtopology(topology, solute_pdb.positions, run.ligand_resname, invert=True)
    ligand_top = _make_subtopology(topology, solute_pdb.positions, run.ligand_resname, invert=False)

    systems = {
        "complex_gb_polar": _build_gb_system(topology, forcefield, include_sa=False),
        "receptor_gb_polar": _build_gb_system(receptor_top, forcefield, include_sa=False),
        "ligand_gb_polar": _build_gb_system(ligand_top, forcefield, include_sa=False),
    }
    contexts = {key: _make_context(system) for key, system in systems.items()}
    h_relax_ctx, h_relax_force, h_relax_heavy = _build_h_relax_context(topology, forcefield)

    universe = mda.Universe(str(run.topology_pdb), str(run.trajectory_dcd))
    snap_idx = _select_snapshot_indices(
        n_frames=len(universe.trajectory),
        discard_fraction=discard_fraction,
        n_snapshots=n_snapshots,
        dt_ps=getattr(universe.trajectory, "dt", None),
        sample_window_ns=sample_window_ns,
    )

    # The analysis PDBs use p66 chain A residue numbers shifted by -3 relative
    # to the manuscript labels: K103->100, Y181->178, Y188->185.
    k103_any = _residue_atom_indices(topology, residue_number=100, chain_id="A")
    k103_sidechain = _atom_indices_by_name(topology, 100, {"CB", "CG", "CD", "CE", "NZ", "OD1", "ND2"})
    y181_any = _residue_atom_indices(topology, residue_number=178, chain_id="A")
    y181_sidechain = _atom_indices_by_name(topology, 178, {"CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH", "SG"})
    y188_sidechain = _atom_indices_by_name(topology, 185, {"CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"})
    ligand_heavy = [
        atom.index
        for atom in topology.atoms()
        if atom.residue.name == run.ligand_resname and (atom.element is None or atom.element.symbol.upper() != "H")
    ]
    ligand_polar = [
        atom.index
        for atom in topology.atoms()
        if atom.residue.name == run.ligand_resname and atom.element is not None and atom.element.symbol.upper() in {"N", "O"}
    ]

    rows: list[dict[str, float | int | str]] = []
    for frame in snap_idx:
        universe.trajectory[int(frame)]
        solute_nm = np.asarray(universe.atoms.positions, dtype=float) / 10.0
        solute_nm = _apply_h_relax(h_relax_ctx, h_relax_force, h_relax_heavy, solute_nm)
        rec_nm = solute_nm[receptor_idx]
        lig_nm = solute_nm[ligand_idx]

        complex_gb = _energy_of(contexts["complex_gb_polar"], solute_nm, 2)
        receptor_gb = _energy_of(contexts["receptor_gb_polar"], rec_nm, 2)
        ligand_gb = _energy_of(contexts["ligand_gb_polar"], lig_nm, 2)

        rows.append(
            {
                "mutation": run.mutation,
                "safe_label": run.safe_label,
                "replicate": run.replicate,
                "frame": int(frame),
                "gb_complex": complex_gb,
                "gb_receptor": receptor_gb,
                "gb_ligand": ligand_gb,
                "binding_dg_gb": complex_gb - receptor_gb - ligand_gb,
                "k103_ligand_min_nm": _min_distance_nm(solute_nm, k103_any, ligand_heavy),
                "k103_sidechain_ligand_polar_min_nm": _min_distance_nm(solute_nm, k103_sidechain, ligand_polar),
                "y181_ligand_min_nm": _min_distance_nm(solute_nm, y181_any, ligand_heavy),
                "y181_sidechain_ligand_min_nm": _min_distance_nm(solute_nm, y181_sidechain, ligand_heavy),
                "y188_sidechain_ligand_min_nm": _min_distance_nm(solute_nm, y188_sidechain, ligand_heavy),
            }
        )

    del h_relax_ctx
    for context in contexts.values():
        del context

    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path("results/analysis/binding_energy/gb_diagnostics"))
    parser.add_argument("--n-snapshots", type=int, default=100)
    parser.add_argument("--discard-fraction", type=float, default=0.25)
    parser.add_argument("--sample-window-ns", type=float, default=None)
    parser.add_argument(
        "--targets",
        type=str,
        default="WT:1,WT:2,WT:3,K103N:1,K103N:2,K103N:3,Y181C:1,Y181C:2,Y181C:3",
        help="Comma-separated mutation:replicate pairs.",
    )
    args = parser.parse_args()

    root = Path.cwd()
    targets = []
    for token in str(args.targets).split(","):
        if not token.strip():
            continue
        mutation, replicate = token.strip().split(":", 1)
        targets.append((mutation.strip(), int(replicate)))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    frames = [
        _evaluate_run(
            _run_input(root, mutation, replicate),
            n_snapshots=args.n_snapshots,
            discard_fraction=args.discard_fraction,
            sample_window_ns=args.sample_window_ns,
        )
        for mutation, replicate in targets
    ]
    snapshot_df = pd.concat(frames, ignore_index=True)
    snapshot_df.to_csv(args.output_dir / "gb_snapshot_terms.csv", index=False)

    summary = (
        snapshot_df.groupby(["mutation", "replicate"], as_index=False)
        .agg(
            gb_complex_mean=("gb_complex", "mean"),
            gb_receptor_mean=("gb_receptor", "mean"),
            gb_ligand_mean=("gb_ligand", "mean"),
            binding_dg_gb_mean=("binding_dg_gb", "mean"),
            binding_dg_gb_std=("binding_dg_gb", "std"),
            k103_ligand_min_nm=("k103_ligand_min_nm", "mean"),
            k103_sidechain_ligand_polar_min_nm=("k103_sidechain_ligand_polar_min_nm", "mean"),
            y181_ligand_min_nm=("y181_ligand_min_nm", "mean"),
            y181_sidechain_ligand_min_nm=("y181_sidechain_ligand_min_nm", "mean"),
            y188_sidechain_ligand_min_nm=("y188_sidechain_ligand_min_nm", "mean"),
        )
        .sort_values(["mutation", "replicate"])
    )
    wt_lookup = summary[summary["mutation"] == "WT"].set_index("replicate")
    rows = []
    for _, row in summary[summary["mutation"] != "WT"].iterrows():
        wt = wt_lookup.loc[int(row["replicate"])]
        out = row.to_dict()
        for column in [
            "gb_complex_mean",
            "gb_receptor_mean",
            "gb_ligand_mean",
            "binding_dg_gb_mean",
            "k103_ligand_min_nm",
            "k103_sidechain_ligand_polar_min_nm",
            "y181_ligand_min_nm",
            "y181_sidechain_ligand_min_nm",
            "y188_sidechain_ligand_min_nm",
        ]:
            out[f"wt_{column}"] = float(wt[column])
            out[f"delta_{column}"] = float(row[column] - wt[column])
        rows.append(out)
    pd.DataFrame(rows).to_csv(args.output_dir / "gb_wt_referenced_summary.csv", index=False)
    summary.to_csv(args.output_dir / "gb_run_summary.csv", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
