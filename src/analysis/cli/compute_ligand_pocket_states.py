#!/usr/bin/env python3
"""Build shared ligand-pocket state assignments across aligned DOR trajectories.

Outputs are written to:
  results/analysis/ligand_pocket_states/{tables,plots,config,medoids}/

This pass focuses on shared, interpretable structural features:
  - pocket C-alpha RMSD to 4NCG over canonical NNIBP residues
  - ligand center of geometry in the aligned pocket frame
  - RMSD of key ligand atoms used in crystal contact definitions
  - crystal-derived protein-ligand contact distances

The intent is to test whether DRMs redistribute occupancy across distinct
ligand-pocket substates rather than merely shifting a single pooled energy.
"""
from __future__ import annotations

import argparse
import json
import logging
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler

from ..pbc import load_mdtraj_trajectory
from ..result_collector import collect_md_results

NNIBP_P66_CANONICAL = [100, 101, 103, 106, 107, 108, 179, 181, 188, 189, 190, 227, 229, 234, 318]
NNIBP_PALM_CANONICAL = [100, 101, 103, 106, 107, 108, 179, 181, 188, 189, 190]
NNIBP_ENTRANCE_CANONICAL = [227, 229, 234, 318]


@dataclass(frozen=True)
class ContactSpec:
    contact_id: str
    category: str
    protein_resid_auth: int
    protein_atom: str
    ligand_atom: str


def _parse_csv_list(text: str) -> list[str]:
    return [tok.strip() for tok in str(text).replace(";", ",").split(",") if tok.strip()]


def _largest_protein_chain_index(topology) -> int:
    best_idx = -1
    best_count = -1
    for chain in topology.chains:
        count = sum(1 for residue in chain.residues if residue.is_protein)
        if count > best_count:
            best_count = count
            best_idx = int(chain.index)
    if best_idx < 0:
        raise ValueError("No protein chain found.")
    return best_idx


def _infer_element(atom_name: str) -> str:
    name = str(atom_name).strip().upper()
    if name.startswith("CL"):
        return "CL"
    if name.startswith("BR"):
        return "BR"
    return name[:1] if name else ""


def _protein_atom_index(topology, chain_idx: int, resseq: int, atom_name: str) -> int | None:
    matches: list[int] = []
    for atom in topology.atoms:
        residue = atom.residue
        chain = residue.chain
        if int(chain.index) != int(chain_idx):
            continue
        if not residue.is_protein:
            continue
        if int(residue.resSeq) != int(resseq):
            continue
        if str(atom.name).strip() != str(atom_name).strip():
            continue
        matches.append(int(atom.index))
    if len(matches) == 1:
        return matches[0]
    return None


def _protein_ca_indices(topology, chain_idx: int, resseq_list: list[int]) -> np.ndarray:
    out: list[int] = []
    for resseq in resseq_list:
        atom_idx = _protein_atom_index(topology, chain_idx, int(resseq), "CA")
        if atom_idx is None:
            raise ValueError(f"Could not map CA for residue {resseq} on chain index {chain_idx}")
        out.append(int(atom_idx))
    return np.asarray(out, dtype=int)


def _protein_residue_index(topology, chain_idx: int, resseq: int) -> int | None:
    matches = [
        int(residue.index)
        for residue in topology.residues
        if residue.is_protein and int(residue.chain.index) == int(chain_idx) and int(residue.resSeq) == int(resseq)
    ]
    if len(matches) == 1:
        return matches[0]
    return None


def _ligand_residue_index(topology, ligand_resname: str) -> int | None:
    matches = [
        int(residue.index)
        for residue in topology.residues
        if str(residue.name).strip() == str(ligand_resname).strip()
    ]
    if len(matches) == 1:
        return matches[0]
    return None


def _find_reference_ligand_atom_index(topology, ligand_resname: str, atom_name: str) -> int | None:
    matches = [
        int(atom.index)
        for atom in topology.atoms
        if str(atom.residue.name).strip() == str(ligand_resname).strip() and str(atom.name).strip() == str(atom_name).strip()
    ]
    if len(matches) == 1:
        return matches[0]
    return None


def _resolve_ligand_atom_index(
    topology,
    xyz_frame0: np.ndarray,
    ligand_resname: str,
    atom_name: str,
    *,
    ref_atom_xyz: np.ndarray | None = None,
    anchor_xyz: np.ndarray | None = None,
) -> int | None:
    ligand_atoms = [atom for atom in topology.atoms if str(atom.residue.name).strip() == str(ligand_resname).strip()]
    if not ligand_atoms:
        return None

    exact = [int(atom.index) for atom in ligand_atoms if str(atom.name).strip() == str(atom_name).strip()]
    if len(exact) == 1:
        return exact[0]

    prefix = [int(atom.index) for atom in ligand_atoms if str(atom.name).strip().startswith(str(atom_name).strip())]
    if len(prefix) == 1:
        return prefix[0]

    candidates = prefix if prefix else []
    if not candidates:
        target_el = _infer_element(atom_name)
        candidates = [
            int(atom.index)
            for atom in ligand_atoms
            if _infer_element(atom.name) == target_el
        ]
    if not candidates:
        return None
    if len(candidates) == 1:
        return candidates[0]

    if ref_atom_xyz is not None:
        d = np.linalg.norm(xyz_frame0[candidates] - np.asarray(ref_atom_xyz, dtype=float), axis=1)
        return int(candidates[int(np.argmin(d))])
    if anchor_xyz is not None:
        d = np.linalg.norm(xyz_frame0[candidates] - np.asarray(anchor_xyz, dtype=float), axis=1)
        return int(candidates[int(np.argmin(d))])
    return int(candidates[0])


def _build_time_axis_ns(traj, total_ns: float | None, frame_indices: np.ndarray) -> np.ndarray:
    if total_ns is not None and np.isfinite(total_ns) and float(total_ns) > 0.0 and traj.n_frames > 1:
        max_frame = max(1, int(traj.n_frames - 1))
        return (frame_indices.astype(float) / float(max_frame)) * float(total_ns)
    dt_ps = getattr(traj, "timestep", None)
    if dt_ps is None:
        dt_ps = getattr(traj, "dt", None)
    if dt_ps is not None and np.isfinite(dt_ps) and float(dt_ps) > 0.0:
        dt_ps = float(dt_ps)
        if dt_ps > 1000.0:
            candidate = dt_ps / 1000.0
            if candidate <= 1000.0:
                dt_ps = candidate
        return frame_indices.astype(float) * dt_ps / 1000.0
    return frame_indices.astype(float)


def _infer_total_ns(row: pd.Series) -> float | None:
    rep_dir = Path(str(row["analysis_dcd"])).parent
    safe = str(row["safe_label"])
    rep = int(row["replicate"])
    state_csv = rep_dir / f"{safe}_rep{rep:02d}_md_state.csv"
    if not state_csv.exists():
        return None
    try:
        df = pd.read_csv(state_csv)
    except Exception:
        return None
    step_col = None
    for candidate in ('#"Step"', "Step"):
        if candidate in df.columns:
            step_col = candidate
            break
    if step_col is None or df.empty:
        return None
    steps = pd.to_numeric(df[step_col], errors="coerce").dropna()
    if steps.empty:
        return None
    return float(steps.max()) * 2.0 / 1_000_000.0


def _aligned_dcd_path(analysis_dcd: Path, aligned_suffix: str) -> Path:
    return analysis_dcd.with_name(f"{analysis_dcd.stem}{aligned_suffix}.dcd")


def _parse_contact_specs(contact_defs: pd.DataFrame) -> list[ContactSpec]:
    return [
        ContactSpec(
            contact_id=str(row["contact_id"]),
            category=str(row["category"]),
            protein_resid_auth=int(row["protein_resid_auth"]),
            protein_atom=str(row["protein_atom"]),
            ligand_atom=str(row["ligand_atom"]),
        )
        for _, row in contact_defs.iterrows()
    ]


def _reference_contact_residue_specs(reference, chain_idx: int, ligand_resname: str, cutoff_angstrom: float) -> list[tuple[int, str]]:
    import mdtraj as md

    ligand_res_idx = _ligand_residue_index(reference.topology, ligand_resname)
    if ligand_res_idx is None:
        raise ValueError(f"Could not locate ligand residue {ligand_resname} in reference structure")

    protein_residues = [
        residue
        for residue in reference.topology.residues
        if residue.is_protein and int(residue.chain.index) == int(chain_idx)
    ]
    contact_pairs = np.asarray([[int(residue.index), int(ligand_res_idx)] for residue in protein_residues], dtype=int)
    distances_nm, returned_pairs = md.compute_contacts(reference, contacts=contact_pairs, scheme="closest-heavy")

    out: list[tuple[int, str]] = []
    for (residue_idx, _lig_idx), distance_nm in zip(returned_pairs.tolist(), distances_nm[0].tolist()):
        residue = reference.topology.residue(int(residue_idx))
        distance_angstrom = float(distance_nm) * 10.0
        if distance_angstrom <= float(cutoff_angstrom):
            label = f"residue_min_distance_{str(residue.name).strip()}{int(residue.resSeq)}_angstrom"
            out.append((int(residue.resSeq), label))
    if not out:
        raise ValueError(f"No crystal-neighbor residues found within {cutoff_angstrom:.2f} A")
    return out


def _sample_frame_indices(n_frames: int, frame_stride: int, max_frames: int | None) -> np.ndarray:
    idx = np.arange(0, int(n_frames), max(1, int(frame_stride)), dtype=int)
    if idx.size == 0 and n_frames > 0:
        idx = np.array([0], dtype=int)
    if idx.size > 1 and idx[-1] != (n_frames - 1):
        idx = np.concatenate([idx, np.array([n_frames - 1], dtype=int)])
    if max_frames is not None and int(max_frames) > 0 and idx.size > int(max_frames):
        take = np.linspace(0, idx.size - 1, num=int(max_frames), dtype=int)
        idx = idx[take]
    return np.unique(idx)


def _plot_state_occupancy_heatmap(occupancy_df: pd.DataFrame, output_png: Path) -> None:
    if occupancy_df.empty:
        return
    pivot = occupancy_df.pivot(index="mutation", columns="state_id", values="occupancy_mean").fillna(0.0)
    muts = pivot.index.tolist()
    states = [int(c) for c in pivot.columns.tolist()]

    fig_w = max(7.5, 0.65 * len(states) + 4.0)
    fig_h = max(5.5, 0.33 * len(muts) + 2.4)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(pivot.to_numpy(dtype=float), aspect="auto", cmap="viridis", vmin=0.0, vmax=max(0.25, float(pivot.to_numpy(dtype=float).max())))
    ax.set_xticks(np.arange(len(states)))
    ax.set_xticklabels([f"State {s}" for s in states], rotation=45, ha="right")
    ax.set_yticks(np.arange(len(muts)))
    ax.set_yticklabels(muts)
    ax.set_title("Mean State Occupancy by Mutation")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cbar.set_label("Occupancy")
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _extract_state_features(
    row: pd.Series,
    *,
    reference,
    ref_pocket_ca_idx: np.ndarray,
    ref_pose_atom_indices: list[int],
    ref_pose_atom_xyz: dict[str, np.ndarray],
    contact_residue_specs: list[tuple[int, str]],
    aligned_suffix: str,
    ligand_resname: str,
    resid_offset: int,
    pocket_resseq: list[int],
    palm_resseq: list[int],
    entrance_resseq: list[int],
    frame_stride: int,
    max_frames: int | None,
) -> tuple[pd.DataFrame, list[str]]:
    analysis_dcd = Path(str(row["analysis_dcd"]))
    topo_path = Path(str(row["analysis_topology_pdb"]))
    aligned_dcd = _aligned_dcd_path(analysis_dcd, aligned_suffix)
    if not aligned_dcd.exists():
        raise FileNotFoundError(f"Missing aligned DCD: {aligned_dcd}")
    if not topo_path.exists():
        raise FileNotFoundError(f"Missing topology: {topo_path}")

    traj = load_mdtraj_trajectory(aligned_dcd, topo_path)
    if traj.n_frames < 1:
        raise ValueError(f"Empty aligned trajectory: {aligned_dcd}")

    frame_idx = _sample_frame_indices(traj.n_frames, frame_stride=frame_stride, max_frames=max_frames)
    sub = traj.slice(frame_idx, copy=True)
    frame0_xyz = traj.xyz[0]

    traj_chain_idx = _largest_protein_chain_index(traj.topology)
    traj_pocket_ca_idx = _protein_ca_indices(
        traj.topology,
        traj_chain_idx,
        [int(auth_resseq) + int(resid_offset) for auth_resseq in pocket_resseq],
    )
    traj_palm_ca_idx = _protein_ca_indices(
        traj.topology,
        traj_chain_idx,
        [int(auth_resseq) + int(resid_offset) for auth_resseq in palm_resseq],
    )
    traj_entrance_ca_idx = _protein_ca_indices(
        traj.topology,
        traj_chain_idx,
        [int(auth_resseq) + int(resid_offset) for auth_resseq in entrance_resseq],
    )

    pose_atom_names = list(ref_pose_atom_xyz.keys())
    traj_pose_atom_idx: list[int] = []
    for atom_name in pose_atom_names:
        atom_idx = _resolve_ligand_atom_index(
            traj.topology,
            xyz_frame0=frame0_xyz,
            ligand_resname=ligand_resname,
            atom_name=atom_name,
            ref_atom_xyz=ref_pose_atom_xyz[atom_name],
        )
        if atom_idx is None:
            raise ValueError(f"Could not map ligand pose atom {atom_name} in {aligned_dcd}")
        traj_pose_atom_idx.append(int(atom_idx))

    import mdtraj as md

    ligand_res_idx = _ligand_residue_index(traj.topology, ligand_resname)
    if ligand_res_idx is None:
        raise ValueError(f"Could not locate ligand residue {ligand_resname} in {aligned_dcd}")

    residue_labels: list[str] = []
    residue_pairs: list[list[int]] = []
    residue_value_map: dict[str, np.ndarray] = {}
    for auth_resseq, label in contact_residue_specs:
        protein_residue_idx = _protein_residue_index(
            traj.topology,
            traj_chain_idx,
            int(auth_resseq) + int(resid_offset),
        )
        if protein_residue_idx is None:
            residue_value_map[label] = np.full(sub.n_frames, np.nan, dtype=float)
            continue
        residue_labels.append(label)
        residue_pairs.append([int(protein_residue_idx), int(ligand_res_idx)])
    if residue_pairs:
        residue_dist_nm, _returned_pairs = md.compute_contacts(
            sub,
            contacts=np.asarray(residue_pairs, dtype=int),
            scheme="closest-heavy",
        )
        residue_dist_angstrom = residue_dist_nm * 10.0
        for idx, label in enumerate(residue_labels):
            residue_value_map[label] = residue_dist_angstrom[:, idx].astype(float)
    if not residue_value_map:
        raise ValueError(f"No residue minimum-distance features could be constructed in {aligned_dcd}")

    ligand_xyz = sub.xyz[:, np.asarray(traj_pose_atom_idx, dtype=int), :]
    ligand_center = ligand_xyz.mean(axis=1)
    pocket_center = sub.xyz[:, traj_pocket_ca_idx, :].mean(axis=1)
    palm_center = sub.xyz[:, traj_palm_ca_idx, :].mean(axis=1)
    entrance_center = sub.xyz[:, traj_entrance_ca_idx, :].mean(axis=1)
    depth_axis = palm_center - entrance_center
    depth_axis_norm = np.linalg.norm(depth_axis, axis=1, keepdims=True)
    depth_axis_unit = depth_axis / np.clip(depth_axis_norm, 1.0e-8, None)
    ligand_from_entrance = ligand_center - entrance_center
    palm_distance = np.linalg.norm(ligand_center - palm_center, axis=1) * 10.0
    entrance_distance = np.linalg.norm(ligand_center - entrance_center, axis=1) * 10.0
    pocket_center_distance = np.linalg.norm(ligand_center - pocket_center, axis=1) * 10.0
    palm_depth_projection = np.sum(ligand_from_entrance * depth_axis_unit, axis=1) * 10.0
    ligand_pose_rmsd = md.rmsd(
        sub,
        reference,
        frame=0,
        atom_indices=np.asarray(traj_pose_atom_idx, dtype=int),
        ref_atom_indices=np.asarray(ref_pose_atom_indices, dtype=int),
    ) * 10.0
    pocket_rmsd = md.rmsd(
        sub,
        reference,
        frame=0,
        atom_indices=np.asarray(traj_pocket_ca_idx, dtype=int),
        ref_atom_indices=np.asarray(ref_pocket_ca_idx, dtype=int),
    ) * 10.0

    total_ns = _infer_total_ns(row)
    time_ns = _build_time_axis_ns(traj, total_ns, frame_idx)

    data: dict[str, object] = {
        "structure": str(row["structure"]),
        "mutation": str(row["mutation"]),
        "safe_label": str(row["safe_label"]),
        "replicate": int(row["replicate"]),
        "frame_index": frame_idx.astype(int),
        "time_ns": time_ns.astype(float),
        "aligned_dcd": str(aligned_dcd),
        "analysis_topology_pdb": str(topo_path),
        "fold_reduction": row.get("fold_reduction", np.nan),
        "pocket_ca_rmsd_angstrom": pocket_rmsd.astype(float),
        "ligand_pose_rmsd_angstrom": ligand_pose_rmsd.astype(float),
        "ligand_palm_distance_angstrom": palm_distance.astype(float),
        "ligand_entrance_distance_angstrom": entrance_distance.astype(float),
        "ligand_pocket_center_distance_angstrom": pocket_center_distance.astype(float),
        "ligand_palm_depth_projection_angstrom": palm_depth_projection.astype(float),
    }
    for _auth_resseq, label in contact_residue_specs:
        data[label] = residue_value_map.get(label, np.full(sub.n_frames, np.nan, dtype=float))

    return pd.DataFrame(data), [label for _auth_resseq, label in contact_residue_specs]


def _save_medoid_frames(frame_df: pd.DataFrame, medoid_df: pd.DataFrame, output_dir: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    output_dir.mkdir(parents=True, exist_ok=True)

    for _, medoid in medoid_df.iterrows():
        dcd_path = Path(str(medoid["aligned_dcd"]))
        topo_path = Path(str(medoid["analysis_topology_pdb"]))
        frame_index = int(medoid["frame_index"])
        state_id = int(medoid["state_id"])
        mutation = str(medoid["mutation"])
        replicate = int(medoid["replicate"])

        traj = load_mdtraj_trajectory(dcd_path, topo_path)
        frame = traj.slice([frame_index], copy=True)
        out_pdb = output_dir / f"state_{state_id:02d}_{mutation.replace('+', '_')}_rep{replicate:02d}_frame{frame_index:04d}.pdb"
        frame.save_pdb(str(out_pdb))
        row = dict(medoid)
        row["medoid_pdb"] = str(out_pdb)
        rows.append(row)

    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Compute shared ligand-pocket states across aligned DOR trajectories.")
    parser.add_argument("--manifest", type=Path, default=Path("manifests/md_manifest.csv"))
    parser.add_argument("--results-dir", type=Path, default=Path("results"))
    parser.add_argument("--reference-cif", type=Path, default=Path("data/structures/4NCG.cif"))
    parser.add_argument(
        "--contact-defs",
        type=Path,
        default=Path("results/tables/holo/dor_key_contact_definitions_4ncg.csv"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/analysis/ligand_pocket_states"),
    )
    parser.add_argument("--aligned-suffix", type=str, default="_aligned_4ncg_ca")
    parser.add_argument("--ligand-resname", type=str, default="2KW")
    parser.add_argument("--resid-offset", type=int, default=-3)
    parser.add_argument("--frame-stride", type=int, default=5)
    parser.add_argument("--max-frames-per-replicate", type=int, default=0)
    parser.add_argument("--n-states", type=int, default=6)
    parser.add_argument("--crystal-contact-cutoff-angstrom", type=float, default=4.0)
    parser.add_argument(
        "--canonical-pocket-resseq",
        type=str,
        default=",".join(str(x) for x in NNIBP_P66_CANONICAL),
    )
    parser.add_argument(
        "--palm-resseq",
        type=str,
        default=",".join(str(x) for x in NNIBP_PALM_CANONICAL),
    )
    parser.add_argument(
        "--entrance-resseq",
        type=str,
        default=",".join(str(x) for x in NNIBP_ENTRANCE_CANONICAL),
    )
    parser.add_argument(
        "--mutations",
        type=str,
        default="",
        help="Optional comma-separated mutation subset.",
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    if not args.reference_cif.exists():
        raise FileNotFoundError(args.reference_cif)
    if not args.contact_defs.exists():
        raise FileNotFoundError(args.contact_defs)

    import mdtraj as md

    run_df = collect_md_results(args.manifest, args.results_dir)
    if run_df.empty:
        raise ValueError("No MD runs found.")

    wanted = set(_parse_csv_list(args.mutations))
    if wanted:
        run_df = run_df[run_df["mutation"].astype(str).isin(wanted)].copy()
        if run_df.empty:
            raise ValueError("Requested mutations did not match any MD runs.")

    pocket_resseq = [int(x) for x in _parse_csv_list(args.canonical_pocket_resseq)]
    palm_resseq = [int(x) for x in _parse_csv_list(args.palm_resseq)]
    entrance_resseq = [int(x) for x in _parse_csv_list(args.entrance_resseq)]
    max_frames = int(args.max_frames_per_replicate) if int(args.max_frames_per_replicate) > 0 else None
    n_states = max(2, int(args.n_states))

    out_tables = args.output_dir / "tables"
    out_plots = args.output_dir / "plots"
    out_config = args.output_dir / "config"
    out_medoids = args.output_dir / "medoids"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_plots.mkdir(parents=True, exist_ok=True)
    out_config.mkdir(parents=True, exist_ok=True)
    out_medoids.mkdir(parents=True, exist_ok=True)

    reference = md.load(str(args.reference_cif))
    ref_chain_idx = _largest_protein_chain_index(reference.topology)
    ref_pocket_ca_idx_arr = _protein_ca_indices(reference.topology, ref_chain_idx, pocket_resseq)

    contact_defs = pd.read_csv(args.contact_defs)
    contact_residue_specs = _reference_contact_residue_specs(
        reference,
        ref_chain_idx,
        args.ligand_resname,
        cutoff_angstrom=float(args.crystal_contact_cutoff_angstrom),
    )

    pose_atom_names = sorted({str(x) for x in contact_defs["ligand_atom"].dropna().astype(str)})
    ref_pose_atom_indices: list[int] = []
    ref_pose_atom_xyz: dict[str, np.ndarray] = {}
    for atom_name in pose_atom_names:
        idx = _find_reference_ligand_atom_index(reference.topology, args.ligand_resname, atom_name)
        if idx is None:
            raise ValueError(f"Could not map reference ligand atom {atom_name} in {args.reference_cif}")
        ref_pose_atom_indices.append(int(idx))
        ref_pose_atom_xyz[str(atom_name)] = reference.xyz[0, int(idx)].copy()

    all_rows: list[pd.DataFrame] = []
    feature_cols: list[str] | None = None
    failures: list[dict[str, object]] = []

    for _, row in run_df.sort_values(["mutation", "replicate"]).iterrows():
        try:
            frame_df, contact_feature_cols = _extract_state_features(
                row,
                reference=reference,
                ref_pocket_ca_idx=ref_pocket_ca_idx_arr,
                ref_pose_atom_indices=ref_pose_atom_indices,
                ref_pose_atom_xyz=ref_pose_atom_xyz,
                contact_residue_specs=contact_residue_specs,
                aligned_suffix=args.aligned_suffix,
                ligand_resname=args.ligand_resname,
                resid_offset=args.resid_offset,
                pocket_resseq=pocket_resseq,
                palm_resseq=palm_resseq,
                entrance_resseq=entrance_resseq,
                frame_stride=args.frame_stride,
                max_frames=max_frames,
            )
            all_rows.append(frame_df)
            if feature_cols is None:
                feature_cols = [
                    "pocket_ca_rmsd_angstrom",
                    "ligand_pose_rmsd_angstrom",
                    "ligand_palm_distance_angstrom",
                    "ligand_entrance_distance_angstrom",
                    "ligand_pocket_center_distance_angstrom",
                    "ligand_palm_depth_projection_angstrom",
                ] + contact_feature_cols
            logging.info(
                "features %s rep%d frames=%d aligned=%s",
                row["mutation"],
                int(row["replicate"]),
                len(frame_df),
                Path(str(frame_df.iloc[0]["aligned_dcd"])).name,
            )
        except Exception as exc:
            failures.append(
                {
                    "mutation": str(row["mutation"]),
                    "replicate": int(row["replicate"]),
                    "error": str(exc),
                }
            )
            logging.warning("state features failed for %s rep%d: %s", row["mutation"], int(row["replicate"]), exc)

    if not all_rows or feature_cols is None:
        raise ValueError("No frame-level features were extracted.")

    frame_df = pd.concat(all_rows, ignore_index=True)
    feature_matrix = frame_df[feature_cols].to_numpy(dtype=float)
    imputer = SimpleImputer(strategy="median")
    scaler = StandardScaler()
    x_imputed = imputer.fit_transform(feature_matrix)
    x_scaled = scaler.fit_transform(x_imputed)

    kmeans = KMeans(n_clusters=n_states, random_state=0, n_init="auto")
    raw_state = kmeans.fit_predict(x_scaled)
    state_counts = pd.Series(raw_state).value_counts().sort_values(ascending=False)
    remap = {int(old): int(new) for new, old in enumerate(state_counts.index.tolist())}
    frame_df["state_id"] = pd.Series(raw_state).map(remap).astype(int)

    centers = kmeans.cluster_centers_
    d2 = np.sum((x_scaled - centers[raw_state]) ** 2, axis=1)
    frame_df["state_distance_to_centroid"] = d2.astype(float)
    frame_df["state_label"] = frame_df["state_id"].map(lambda s: f"State {int(s)}")

    occupancy_rep_raw = (
        frame_df.groupby(["mutation", "replicate", "state_id"], as_index=False)
        .agg(n_frames=("frame_index", "count"))
    )
    rep_totals = occupancy_rep_raw.groupby(["mutation", "replicate"])["n_frames"].sum().rename("n_frames_total")
    rep_keys = frame_df[["mutation", "replicate"]].drop_duplicates().sort_values(["mutation", "replicate"]).reset_index(drop=True)
    state_keys = pd.DataFrame({"state_id": sorted(frame_df["state_id"].dropna().astype(int).unique().tolist())})
    occupancy_rep = (
        rep_keys.assign(_tmp=1)
        .merge(state_keys.assign(_tmp=1), on="_tmp", how="inner")
        .drop(columns="_tmp")
    )
    occupancy_rep = occupancy_rep.merge(
        occupancy_rep_raw,
        on=["mutation", "replicate", "state_id"],
        how="left",
    )
    occupancy_rep["n_frames"] = occupancy_rep["n_frames"].fillna(0).astype(int)
    occupancy_rep = occupancy_rep.merge(rep_totals, on=["mutation", "replicate"], how="left")
    occupancy_rep["occupancy"] = occupancy_rep["n_frames"] / occupancy_rep["n_frames_total"].replace(0, np.nan)

    occupancy_mut = (
        occupancy_rep.groupby(["mutation", "state_id"], as_index=False)
        .agg(
            occupancy_mean=("occupancy", "mean"),
            occupancy_std=("occupancy", "std"),
            occupancy_sem=("occupancy", lambda x: float(np.std(x, ddof=1) / np.sqrt(len(x))) if len(x) > 1 else 0.0),
            n_replicates=("replicate", "nunique"),
        )
    )

    transitions: list[dict[str, object]] = []
    for (mutation, replicate), grp in frame_df.sort_values(["mutation", "replicate", "frame_index"]).groupby(["mutation", "replicate"]):
        states = grp["state_id"].to_numpy(dtype=int)
        if states.size < 2:
            continue
        for a, b in zip(states[:-1], states[1:]):
            transitions.append(
                {
                    "mutation": str(mutation),
                    "replicate": int(replicate),
                    "state_from": int(a),
                    "state_to": int(b),
                    "count": 1,
                }
            )
    transition_df = (
        pd.DataFrame(transitions)
        .groupby(["mutation", "replicate", "state_from", "state_to"], as_index=False)
        .agg(count=("count", "sum"))
        if transitions
        else pd.DataFrame(columns=["mutation", "replicate", "state_from", "state_to", "count"])
    )

    state_summary = frame_df.groupby("state_id")[feature_cols].agg(["mean", "std"])
    state_summary.columns = [f"{col[0]}_{col[1]}" for col in state_summary.columns.to_flat_index()]
    state_summary = state_summary.reset_index()
    state_counts_df = frame_df.groupby("state_id", as_index=False).size().rename(columns={"size": "n_frames"})
    state_summary = state_summary.merge(state_counts_df, on="state_id", how="left")
    state_summary["occupancy_all_frames"] = state_summary["n_frames"] / float(len(frame_df))

    medoid_rows: list[dict[str, object]] = []
    for state_id, grp in frame_df.groupby("state_id"):
        medoid = grp.sort_values("state_distance_to_centroid").iloc[0]
        medoid_rows.append(dict(medoid))
    medoid_df = pd.DataFrame(medoid_rows).sort_values("state_id").reset_index(drop=True)
    medoid_df = _save_medoid_frames(frame_df, medoid_df, out_medoids)

    frame_df.to_csv(out_tables / "frame_features.csv", index=False)
    frame_df[
        [
            "mutation",
            "safe_label",
            "replicate",
            "frame_index",
            "time_ns",
            "aligned_dcd",
            "analysis_topology_pdb",
            "state_id",
            "state_label",
            "state_distance_to_centroid",
        ]
    ].to_csv(out_tables / "state_assignments.csv", index=False)
    occupancy_rep.to_csv(out_tables / "state_occupancy_by_replicate.csv", index=False)
    occupancy_mut.to_csv(out_tables / "state_occupancy_by_mutation.csv", index=False)
    transition_df.to_csv(out_tables / "state_transition_counts.csv", index=False)
    state_summary.to_csv(out_tables / "state_summary.csv", index=False)
    medoid_df.to_csv(out_tables / "state_medoids.csv", index=False)
    pd.DataFrame(failures).to_csv(out_tables / "feature_extraction_failures.csv", index=False)

    config = {
        "manifest": str(args.manifest),
        "results_dir": str(args.results_dir),
        "reference_cif": str(args.reference_cif),
        "contact_defs": str(args.contact_defs),
        "crystal_contact_cutoff_angstrom": float(args.crystal_contact_cutoff_angstrom),
        "crystal_contact_residue_specs": [
            {"auth_resseq": int(auth_resseq), "label": str(label)}
            for auth_resseq, label in contact_residue_specs
        ],
        "aligned_suffix": str(args.aligned_suffix),
        "ligand_resname": str(args.ligand_resname),
        "resid_offset": int(args.resid_offset),
        "frame_stride": int(args.frame_stride),
        "max_frames_per_replicate": None if max_frames is None else int(max_frames),
        "n_states": int(n_states),
        "canonical_pocket_resseq": [int(x) for x in pocket_resseq],
        "palm_resseq": [int(x) for x in palm_resseq],
        "entrance_resseq": [int(x) for x in entrance_resseq],
        "feature_columns": list(feature_cols),
        "n_replicates_succeeded": int(frame_df[["mutation", "replicate"]].drop_duplicates().shape[0]),
        "n_replicates_failed": int(len(failures)),
        "n_frames_total": int(len(frame_df)),
    }
    (out_config / "run_config.json").write_text(json.dumps(config, indent=2))

    _plot_state_occupancy_heatmap(occupancy_mut, out_plots / "state_occupancy_by_mutation.png")

    logging.info(
        "state analysis complete: frames=%d replicates=%d states=%d failures=%d output=%s",
        len(frame_df),
        frame_df[["mutation", "replicate"]].drop_duplicates().shape[0],
        n_states,
        len(failures),
        args.output_dir,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
