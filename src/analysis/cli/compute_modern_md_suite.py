#!/usr/bin/env python3
"""Modern NNIBP descriptor suite from pre-imaged analysis DCDs.

Computes (standard libraries only — thin I/O wrappers):
  - DOR–protein H-bond occupancy by residue   (MDAnalysis HydrogenBondAnalysis)
  - Ligand RMSF after NNIBP Cα fit            (mdtraj.rmsf)
  - Pose clusters on aligned DOR heavy atoms  (scipy.cluster.vq.kmeans2)
  - NNIBP pocket volume proxy                 (src.analysis.metrics)
  - NNIBP Cα PCA across the panel             (numpy SVD)
  - NNIBP–NNIBP DCCM + WT-difference heatmaps (Pearson on Cα)
  - NNIBP residue–residue contact networks    (heavy-atom cutoff frequency)

Requires *_analysis_pbcfix.dcd (write once with fix_pbc_trajectories).

    ~/miniconda3/envs/nnrti-prep/bin/python -m src.analysis.cli.compute_modern_md_suite
    ~/miniconda3/envs/nnrti-prep/bin/python -m src.analysis.cli.compute_modern_md_suite --from-tables
"""
from __future__ import annotations

import argparse
import json
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..md_timing import infer_production_ns
from ..metrics import pocket_volume_proxy_from_universe
from ..pbc import load_mdtraj_trajectory, pbcfix_dcd_for, raw_analysis_dcd_for
from ..result_collector import _prepare_profile_jobs, collect_md_results

REPO = Path(__file__).resolve().parents[3]
LOGGER = logging.getLogger("modern_md_suite")

RESID_OFFSET = -3
NNIBP_AUTH = (100, 101, 103, 106, 107, 108, 179, 181, 188, 189, 190, 227, 229, 234, 318)
STORY_MUTATIONS = ("V106A", "V106I", "V106I+F227C", "V106A+F227L", "G190E", "Y188L")
CONTACT_CUTOFF_A = 4.0
HBOND_DA_CUTOFF = 3.5
HBOND_ANGLE = 135.0


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compute modern NNIBP MD descriptor suite.")
    p.add_argument("--manifest", type=Path, default=Path("manifests/md_manifest.csv"))
    p.add_argument("--output-dir", type=Path, default=Path("results/analysis/modern_md_suite"))
    p.add_argument("--resid-offset", type=int, default=RESID_OFFSET)
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--frame-stride", type=int, default=1, help="Stride for RMSF/PCA/DCCM/contacts")
    p.add_argument("--pocket-max-frames", type=int, default=40, help="Subsample for pocket volume")
    p.add_argument("--hbond-max-frames", type=int, default=80, help="Subsample for H-bond analysis")
    p.add_argument("--n-pose-clusters", type=int, default=4)
    p.add_argument("--mutations", nargs="*", default=None)
    p.add_argument("--from-tables", action="store_true", help="Replot/aggregate from saved tables/npy")
    return p.parse_args()


def _display_mutation(raw: str) -> str:
    text = str(raw).strip()
    return "WT" if text.lower() == "wt" else text


def _mutation_sort_key(mutation: str) -> tuple[int, str]:
    if mutation == "WT":
        return (0, mutation)
    if "+" in mutation:
        return (2, mutation)
    return (1, mutation)


def _nnibp_traj_resseqs(resid_offset: int) -> list[int]:
    return [int(a) + int(resid_offset) for a in NNIBP_AUTH]


def _select_nnibp_ca(topology, resid_offset: int) -> np.ndarray:
    resseqs = _nnibp_traj_resseqs(resid_offset)
    query = " or ".join(f"resSeq {r}" for r in resseqs)
    idx = np.asarray(topology.select(f"chainid 0 and name CA and ({query})"), dtype=int)
    if idx.size != len(NNIBP_AUTH):
        raise ValueError(f"NNIBP Cα count {idx.size} != {len(NNIBP_AUTH)}")
    # Order to match NNIBP_AUTH
    by_res: dict[int, int] = {}
    for atom in topology.atoms:
        if int(atom.index) in set(idx.tolist()):
            by_res[int(atom.residue.resSeq)] = int(atom.index)
    ordered = [by_res[r] for r in resseqs]
    return np.asarray(ordered, dtype=int)


def _frame_subset(n_frames: int, max_frames: int, stride: int) -> np.ndarray:
    base = np.arange(0, n_frames, max(1, stride), dtype=int)
    if max_frames is not None and len(base) > max_frames:
        pick = np.linspace(0, len(base) - 1, num=max_frames, dtype=int)
        base = base[pick]
    return base


def _process_job(
    job: dict,
    resid_offset: int,
    frame_stride: int,
    pocket_max_frames: int,
    hbond_max_frames: int,
    contact_cutoff: float,
    npy_dir: str,
) -> dict:
    """Per-replicate descriptors. Returns serializable dict (+ writes npy sidecars)."""
    import mdtraj as md
    import MDAnalysis as mda
    from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
    from MDAnalysis.lib.distances import capped_distance

    mutation = _display_mutation(job["mutation"])
    replicate = int(job["replicate"])
    production_ps = float(job["production_ps"])
    traj_path = Path(job["trajectory"])
    topo_path = Path(job["topology"])
    traj = load_mdtraj_trajectory(traj_path, topo_path)
    n_frames = int(traj.n_frames)
    if n_frames < 2:
        raise ValueError("need ≥2 frames")

    lig_idx = np.asarray(traj.topology.select("resname '2KW' and not element H"), dtype=int)
    if lig_idx.size == 0:
        raise ValueError("no DOR heavy atoms")
    nnibp_ca = _select_nnibp_ca(traj.topology, resid_offset)

    # --- Ligand RMSF after NNIBP Cα fit to mean ---
    pose = traj[:: max(1, frame_stride)]
    pose.superpose(pose, frame=0, atom_indices=nnibp_ca)
    # Manual RMSF in the pocket-aligned frame (md.rmsf mis-handles this setup).
    lig_xyz_A = pose.xyz[:, lig_idx, :] * 10.0
    lig_rmsf = np.sqrt(((lig_xyz_A - lig_xyz_A.mean(axis=0)) ** 2).sum(axis=-1).mean(axis=0))
    lig_names = [traj.topology.atom(i).name for i in lig_idx]
    rmsf_rows = [
        {
            "mutation": mutation,
            "replicate": replicate,
            "atom_name": name,
            "atom_index": int(idx),
            "rmsf_angstrom": float(val),
        }
        for name, idx, val in zip(lig_names, lig_idx, lig_rmsf)
    ]

    # Save aligned ligand xyz (for pose clustering) and NNIBP CA xyz (PCA/DCCM)
    # Use a common stride subset
    idx_sub = _frame_subset(n_frames, max_frames=120, stride=frame_stride)
    sub = traj[idx_sub]
    sub.superpose(sub, frame=0, atom_indices=nnibp_ca)
    nnibp_xyz = sub.xyz[:, nnibp_ca, :] * 10.0  # Å
    lig_xyz = sub.xyz[:, lig_idx, :] * 10.0
    safe = str(job["safe_label"])
    npy_root = Path(npy_dir)
    npy_root.mkdir(parents=True, exist_ok=True)
    stem = f"{safe}_rep{replicate:02d}"
    np.save(npy_root / f"{stem}_nnibp_ca_xyz.npy", nnibp_xyz.astype(np.float32))
    np.save(npy_root / f"{stem}_lig_heavy_xyz.npy", lig_xyz.astype(np.float32))
    np.save(npy_root / f"{stem}_frame_indices.npy", idx_sub.astype(np.int32))

    # --- NNIBP–NNIBP contact frequency + per-frame DCCM ingredients ---
    # Contact: heavy atoms of each NNIBP residue pair within cutoff
    nnibp_heavy = []
    for auth, t_res in zip(NNIBP_AUTH, _nnibp_traj_resseqs(resid_offset)):
        h = np.asarray(
            traj.topology.select(f"chainid 0 and resSeq {t_res} and not element H"),
            dtype=int,
        )
        if h.size == 0:
            raise ValueError(f"no heavy atoms for auth {auth}")
        nnibp_heavy.append(h)
    n_res = len(NNIBP_AUTH)
    contact_counts = np.zeros((n_res, n_res), dtype=float)
    for fi in idx_sub:
        xyz = traj.xyz[fi] * 10.0
        for i in range(n_res):
            for j in range(i + 1, n_res):
                a = xyz[nnibp_heavy[i]]
                b = xyz[nnibp_heavy[j]]
                d2 = ((a[:, None, :] - b[None, :, :]) ** 2).sum(axis=-1)
                if np.any(d2 < contact_cutoff**2):
                    contact_counts[i, j] += 1.0
                    contact_counts[j, i] += 1.0
    contact_freq = contact_counts / float(len(idx_sub))
    np.save(npy_root / f"{stem}_nnibp_contact_freq.npy", contact_freq.astype(np.float32))

    # --- Pocket volume (subsampled) ---
    pocket_idx = _frame_subset(n_frames, max_frames=pocket_max_frames, stride=max(1, n_frames // pocket_max_frames))
    # metrics.pocket_volume_proxy loads full traj internally — call per-frame via universe
    u = mda.Universe(str(topo_path), str(traj_path))

    pocket_vals = []
    for fi in pocket_idx:
        u.trajectory[int(fi)]
        pocket_vals.append(float(pocket_volume_proxy_from_universe(u)))
    pocket_vals = np.asarray(pocket_vals, dtype=float)

    # --- H-bonds DOR↔protein (evenly strided frames) ---
    hb_step = max(1, n_frames // max(1, hbond_max_frames))
    hb_idx = np.arange(0, n_frames, hb_step, dtype=int)
    hbond = HydrogenBondAnalysis(
        u,
        donors_sel="resname 2KW or protein",
        hydrogens_sel="name H*",
        acceptors_sel="resname 2KW or protein",
        between=["resname 2KW", "protein"],
        d_a_cutoff=HBOND_DA_CUTOFF,
        d_h_a_angle_cutoff=HBOND_ANGLE,
    )
    hbond.run(start=0, stop=n_frames, step=hb_step)
    hb = hbond.results.hbonds

    frame_set = set(int(x) for x in hb_idx.tolist())
    # hbonds columns: frame, donor_idx, hydrogen_idx, acceptor_idx, ...
    residue_frame_pairs: dict[int, set[int]] = {}
    total_per_frame: dict[int, set[tuple[int, int]]] = {f: set() for f in frame_set}
    if hb is not None and len(hb) > 0:
        lig_indices = set(int(i) for i in u.select_atoms("resname 2KW").indices)
        for row in hb:
            frame_id = int(row[0])
            if frame_id not in frame_set:
                continue
            d_idx, a_idx = int(row[1]), int(row[3])
            total_per_frame[frame_id].add((d_idx, a_idx))
            d_atom = u.atoms[d_idx]
            a_atom = u.atoms[a_idx]
            # protein partner residue (traj resid → auth)
            for atom in (d_atom, a_atom):
                if int(atom.index) in lig_indices:
                    continue
                if atom.resname == "2KW":
                    continue
                traj_resid = int(atom.resid)
                auth = traj_resid - int(resid_offset)
                residue_frame_pairs.setdefault(auth, set()).add(frame_id)

    n_hb_frames = float(len(frame_set))
    hbond_rows = []
    for auth, frames in sorted(residue_frame_pairs.items()):
        hbond_rows.append(
            {
                "mutation": mutation,
                "replicate": replicate,
                "auth_resid": int(auth),
                "n_frames_with_hbond": int(len(frames)),
                "occupancy": float(len(frames) / n_hb_frames),
            }
        )
    mean_hbonds = float(np.mean([len(v) for v in total_per_frame.values()])) if total_per_frame else 0.0

    # Ligand COM feature for clustering (already aligned in lig_xyz)
    lig_com = lig_xyz.mean(axis=1)  # (n_sub, 3)

    return {
        "mutation": mutation,
        "safe_label": safe,
        "replicate": replicate,
        "n_frames": n_frames,
        "total_ns": float(production_ps / 1000.0),
        "timing_source": job.get("timing_source"),
        "rmsf_mean_angstrom": float(np.mean(lig_rmsf)),
        "rmsf_max_angstrom": float(np.max(lig_rmsf)),
        "pocket_volume_mean": float(np.mean(pocket_vals)),
        "pocket_volume_std": float(np.std(pocket_vals)),
        "pocket_n_frames": int(len(pocket_vals)),
        "hbond_count_mean": mean_hbonds,
        "hbond_n_frames": int(len(frame_set)),
        "npy_stem": stem,
        "rmsf_rows": rmsf_rows,
        "hbond_rows": hbond_rows,
        "lig_com_mean": lig_com.mean(axis=0).tolist(),
        "lig_com_std": lig_com.std(axis=0).tolist(),
    }


def _dccm_from_xyz(xyz: np.ndarray) -> np.ndarray:
    """xyz: (n_frames, n_atoms, 3) → (n_atoms, n_atoms) Pearson DCCM."""
    n_frames, n_atoms, _ = xyz.shape
    mean = xyz.mean(axis=0)
    disp = xyz - mean
    # Flatten xyz to (n_frames, n_atoms*3) then per-atom vector correlation
    # Standard DCCM: C_ij = <Δr_i · Δr_j> / (sqrt<Δr_i²> sqrt<Δr_j²>)
    flat = disp.reshape(n_frames, n_atoms, 3)
    corr = np.zeros((n_atoms, n_atoms), dtype=float)
    norms = np.sqrt((flat**2).sum(axis=(0, 2)) / n_frames)
    for i in range(n_atoms):
        for j in range(i, n_atoms):
            dot = (flat[:, i, :] * flat[:, j, :]).sum(axis=1).mean()
            denom = norms[i] * norms[j]
            val = float(dot / denom) if denom > 1e-12 else 0.0
            corr[i, j] = corr[j, i] = val
    return corr


def _pca_project(stacked: np.ndarray, n_components: int = 3) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """stacked (n_samples, n_features) → scores, components, explained_var_ratio."""
    x = stacked - stacked.mean(axis=0, keepdims=True)
    # economy SVD
    u, s, vt = np.linalg.svd(x, full_matrices=False)
    components = vt[:n_components]
    scores = u[:, :n_components] * s[:n_components]
    var = (s**2) / max(1, stacked.shape[0] - 1)
    ratio = var / var.sum()
    return scores, components, ratio[:n_components]


def _aggregate_and_plot(out: Path, metas: list[dict], *, n_pose_clusters: int) -> None:
    tables = out / "tables"
    plots = out / "plots"
    npy = out / "npy"
    tables.mkdir(parents=True, exist_ok=True)
    plots.mkdir(parents=True, exist_ok=True)

    inv = pd.DataFrame([{k: m[k] for k in (
        "mutation", "safe_label", "replicate", "n_frames", "total_ns", "timing_source",
        "rmsf_mean_angstrom", "rmsf_max_angstrom", "pocket_volume_mean", "pocket_volume_std",
        "hbond_count_mean", "hbond_n_frames", "npy_stem",
    )} for m in metas])
    inv.to_csv(tables / "replicate_inventory.csv", index=False)

    rmsf_df = pd.DataFrame([r for m in metas for r in m["rmsf_rows"]])
    rmsf_df.to_csv(tables / "ligand_rmsf_per_atom.csv", index=False)
    rmsf_mean = (
        rmsf_df.groupby(["mutation", "atom_name"], as_index=False)["rmsf_angstrom"]
        .agg(rmsf_mean="mean", rmsf_sem=lambda s: float(s.std(ddof=1) / np.sqrt(len(s))) if len(s) > 1 else 0.0)
    )
    rmsf_mean.to_csv(tables / "ligand_rmsf_genotype_mean.csv", index=False)

    hb_delta = None
    hb_df = pd.DataFrame([r for m in metas for r in m["hbond_rows"]])
    if not hb_df.empty:
        hb_df.to_csv(tables / "dor_hbond_occupancy_per_rep.csv", index=False)
        hb_g = (
            hb_df.groupby(["mutation", "auth_resid"], as_index=False)["occupancy"]
            .agg(occupancy_mean="mean", occupancy_sem=lambda s: float(s.std(ddof=1) / np.sqrt(len(s))) if len(s) > 1 else 0.0, n_reps="count")
        )
        hb_g.to_csv(tables / "dor_hbond_occupancy_genotype.csv", index=False)
        wt = hb_g[hb_g["mutation"] == "WT"][["auth_resid", "occupancy_mean"]].rename(
            columns={"occupancy_mean": "wt_occupancy"}
        )
        hb_delta = hb_g.merge(wt, on="auth_resid", how="left")
        hb_delta["delta_vs_wt"] = hb_delta["occupancy_mean"] - hb_delta["wt_occupancy"]
        hb_delta.to_csv(tables / "dor_hbond_occupancy_delta_vs_wt.csv", index=False)

    pocket = inv[["mutation", "replicate", "pocket_volume_mean", "pocket_volume_std", "total_ns"]].copy()
    pocket.to_csv(tables / "pocket_volume_per_rep.csv", index=False)
    pocket_g = (
        pocket.groupby("mutation", as_index=False)["pocket_volume_mean"]
        .agg(volume_mean="mean", volume_sem=lambda s: float(s.std(ddof=1) / np.sqrt(len(s))) if len(s) > 1 else 0.0)
    )
    pocket_g.to_csv(tables / "pocket_volume_genotype.csv", index=False)

    # --- PCA on mean NNIBP CA structure per replicate ---
    feat_rows = []
    labels = []
    for m in metas:
        xyz = np.load(npy / f"{m['npy_stem']}_nnibp_ca_xyz.npy")  # (F, 15, 3)
        mean_struct = xyz.mean(axis=0).ravel()  # 45
        feat_rows.append(mean_struct)
        labels.append((m["mutation"], m["replicate"]))
    feat = np.vstack(feat_rows)
    scores, components, ratio = _pca_project(feat, n_components=3)
    pca_df = pd.DataFrame(
        {
            "mutation": [a for a, _ in labels],
            "replicate": [b for _, b in labels],
            "PC1": scores[:, 0],
            "PC2": scores[:, 1],
            "PC3": scores[:, 2],
        }
    )
    pca_df.to_csv(tables / "nnibp_pca_scores.csv", index=False)
    pd.DataFrame(
        {
            "component": [f"PC{i+1}" for i in range(len(ratio))],
            "explained_variance_ratio": ratio,
        }
    ).to_csv(tables / "nnibp_pca_variance.csv", index=False)
    np.save(npy / "nnibp_pca_components.npy", components)

    # --- DCCM per replicate + WT mean difference for stories ---
    dccm_by_key: dict[tuple[str, int], np.ndarray] = {}
    for m in metas:
        xyz = np.load(npy / f"{m['npy_stem']}_nnibp_ca_xyz.npy")
        dccm = _dccm_from_xyz(xyz)
        dccm_by_key[(m["mutation"], int(m["replicate"]))] = dccm
        np.save(npy / f"{m['npy_stem']}_nnibp_dccm.npy", dccm.astype(np.float32))

    wt_mats = [dccm_by_key[k] for k in dccm_by_key if k[0] == "WT"]
    if wt_mats:
        wt_mean = np.mean(np.stack(wt_mats), axis=0)
        np.save(npy / "wt_mean_nnibp_dccm.npy", wt_mean.astype(np.float32))
        diff_rows = []
        for mut in STORY_MUTATIONS:
            mats = [dccm_by_key[k] for k in dccm_by_key if k[0] == mut]
            if not mats:
                continue
            mut_mean = np.mean(np.stack(mats), axis=0)
            delta = mut_mean - wt_mean
            np.save(npy / f"dccm_delta_{mut.replace('+', '_')}_minus_wt.npy", delta.astype(np.float32))
            for i, ai in enumerate(NNIBP_AUTH):
                for j, aj in enumerate(NNIBP_AUTH):
                    if j < i:
                        continue
                    diff_rows.append(
                        {
                            "mutation": mut,
                            "auth_i": int(ai),
                            "auth_j": int(aj),
                            "dccm_mut": float(mut_mean[i, j]),
                            "dccm_wt": float(wt_mean[i, j]),
                            "delta": float(delta[i, j]),
                        }
                    )
        if diff_rows:
            pd.DataFrame(diff_rows).to_csv(tables / "nnibp_dccm_delta_vs_wt.csv", index=False)

    # --- Contact network Δ vs WT ---
    contact_by_key = {}
    for m in metas:
        contact_by_key[(m["mutation"], int(m["replicate"]))] = np.load(
            npy / f"{m['npy_stem']}_nnibp_contact_freq.npy"
        )
    wt_c = [contact_by_key[k] for k in contact_by_key if k[0] == "WT"]
    net_rows = []
    if wt_c:
        wt_c_mean = np.mean(np.stack(wt_c), axis=0)
        for mut, group in inv.groupby("mutation"):
            mats = [contact_by_key[(mut, int(r))] for r in group["replicate"]]
            mut_mean = np.mean(np.stack(mats), axis=0)
            delta = mut_mean - wt_c_mean
            for i, ai in enumerate(NNIBP_AUTH):
                for j, aj in enumerate(NNIBP_AUTH):
                    if j <= i:
                        continue
                    net_rows.append(
                        {
                            "mutation": mut,
                            "auth_i": int(ai),
                            "auth_j": int(aj),
                            "freq_mut": float(mut_mean[i, j]),
                            "freq_wt": float(wt_c_mean[i, j]),
                            "delta_freq": float(delta[i, j]),
                        }
                    )
        pd.DataFrame(net_rows).to_csv(tables / "nnibp_contact_network_delta_vs_wt.csv", index=False)

    # --- Pose clustering (pooled ligand COM + flattened heavy coords subsample) ---
    from scipy.cluster.vq import kmeans2, whiten

    pose_feats = []
    pose_meta = []
    for m in metas:
        lig = np.load(npy / f"{m['npy_stem']}_lig_heavy_xyz.npy")  # (F, L, 3)
        # per-frame: COM + first 3 principal atom coords relative to COM
        com = lig.mean(axis=1)
        centered = lig - com[:, None, :]
        # take mean abs deviation per atom as compact flexibility fingerprint + COM
        mad = np.mean(np.abs(centered), axis=1)  # (F, L)
        # subsample frames for clustering
        n = lig.shape[0]
        pick = np.linspace(0, n - 1, num=min(40, n), dtype=int)
        for fi in pick:
            feat = np.concatenate([com[fi], mad[fi]])
            pose_feats.append(feat)
            pose_meta.append({"mutation": m["mutation"], "replicate": m["replicate"], "frame_i": int(fi)})
    pose_X = np.vstack(pose_feats)
    pose_Xw = whiten(pose_X)
    centroids, labels_k = kmeans2(pose_Xw, n_pose_clusters, minit="points", seed=0)
    pose_df = pd.DataFrame(pose_meta)
    pose_df["cluster"] = labels_k.astype(int)
    pose_df.to_csv(tables / "pose_cluster_assignments.csv", index=False)
    pop = (
        pose_df.groupby(["mutation", "replicate", "cluster"], as_index=False)
        .size()
        .rename(columns={"size": "n_frames"})
    )
    tot = pop.groupby(["mutation", "replicate"], as_index=False)["n_frames"].sum().rename(columns={"n_frames": "n_total"})
    pop = pop.merge(tot, on=["mutation", "replicate"])
    pop["fraction"] = pop["n_frames"] / pop["n_total"]
    pop.to_csv(tables / "pose_cluster_populations.csv", index=False)
    pop_g = (
        pop.groupby(["mutation", "cluster"], as_index=False)["fraction"]
        .agg(fraction_mean="mean", fraction_sem=lambda s: float(s.std(ddof=1) / np.sqrt(len(s))) if len(s) > 1 else 0.0)
    )
    pop_g.to_csv(tables / "pose_cluster_populations_genotype.csv", index=False)

    # ===== Plots =====
    _plot_rmsf(rmsf_mean, plots / "ligand_rmsf_by_genotype.png")
    _plot_pocket(pocket_g, plots / "pocket_volume_by_genotype.png")
    _plot_pca(pca_df, ratio, plots / "nnibp_pca_pc1_pc2.png")
    if hb_delta is not None:
        _plot_hbond_heatmap(hb_delta, plots / "dor_hbond_delta_heatmap.png")
    _plot_dccm_deltas(npy, plots)
    _plot_pose_pops(pop_g, plots / "pose_cluster_populations.png")
    if net_rows:
        _plot_network_deltas(pd.DataFrame(net_rows), plots / "nnibp_contact_delta_stories.png")


def _plot_rmsf(rmsf_mean: pd.DataFrame, out: Path) -> None:
    focus = [m for m in ("WT", "Y188L", "V106I+F227C", "V106A", "G190E") if m in set(rmsf_mean["mutation"])]
    if not focus:
        focus = sorted(rmsf_mean["mutation"].unique(), key=_mutation_sort_key)[:6]
    fig, ax = plt.subplots(figsize=(10, 4.2))
    atoms = list(dict.fromkeys(rmsf_mean["atom_name"].tolist()))
    x = np.arange(len(atoms))
    for mut in focus:
        sub = rmsf_mean[rmsf_mean["mutation"] == mut].set_index("atom_name").reindex(atoms)
        ax.plot(x, sub["rmsf_mean"], marker="o", ms=3, lw=1.4, label=mut)
    ax.set_xticks(x[:: max(1, len(atoms) // 12)])
    ax.set_xticklabels([atoms[i] for i in range(0, len(atoms), max(1, len(atoms) // 12))], rotation=90, fontsize=7)
    ax.set_ylabel("DOR RMSF (Å)")
    ax.set_title("Ligand RMSF after NNIBP Cα alignment")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(alpha=0.25, linestyle=":")
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def _plot_pocket(pocket_g: pd.DataFrame, out: Path) -> None:
    df = pocket_g.sort_values("mutation", key=lambda s: s.map(lambda m: _mutation_sort_key(m)))
    fig, ax = plt.subplots(figsize=(10, 4))
    x = np.arange(len(df))
    ax.bar(x, df["volume_mean"], yerr=df["volume_sem"], capsize=2, color="#2c6fbb", ecolor="0.4")
    ax.set_xticks(x)
    ax.set_xticklabels(df["mutation"], rotation=90, fontsize=8)
    ax.set_ylabel("NNIBP pocket volume proxy (Å³)")
    ax.set_title("Pocket volume by genotype (mean ± SEM, n = 3)")
    ax.grid(alpha=0.25, axis="y", linestyle=":")
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def _plot_pca(pca_df: pd.DataFrame, ratio: np.ndarray, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    muts = sorted(pca_df["mutation"].unique(), key=_mutation_sort_key)
    cmap = plt.cm.tab20(np.linspace(0, 1, len(muts)))
    for color, mut in zip(cmap, muts):
        sub = pca_df[pca_df["mutation"] == mut]
        ax.scatter(sub["PC1"], sub["PC2"], s=36, color=color, label=mut, alpha=0.85, edgecolors="0.2", linewidths=0.3)
    ax.set_xlabel(f"PC1 ({100*ratio[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({100*ratio[1]:.1f}%)")
    ax.set_title("NNIBP Cα PCA (per-replicate mean structures)")
    ax.legend(fontsize=6, ncol=2, loc="best")
    ax.grid(alpha=0.25, linestyle=":")
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def _plot_hbond_heatmap(hb_delta: pd.DataFrame, out: Path) -> None:
    if "delta_vs_wt" not in hb_delta.columns:
        return
    # Focus on residues that H-bond in WT or any mutant with occupancy > 0.05
    keep_res = sorted(set(hb_delta.loc[hb_delta["occupancy_mean"] > 0.05, "auth_resid"].astype(int)))
    muts = sorted(hb_delta["mutation"].unique(), key=_mutation_sort_key)
    muts = [m for m in muts if m != "WT"]
    if not keep_res or not muts:
        return
    mat = np.full((len(muts), len(keep_res)), np.nan)
    for i, mut in enumerate(muts):
        sub = hb_delta[hb_delta["mutation"] == mut].set_index("auth_resid")
        for j, r in enumerate(keep_res):
            if r in sub.index:
                mat[i, j] = float(sub.loc[r, "delta_vs_wt"])
    fig, ax = plt.subplots(figsize=(max(8, 0.45 * len(keep_res)), max(5, 0.35 * len(muts))))
    lim = float(np.nanmax(np.abs(mat))) if np.isfinite(mat).any() else 1.0
    lim = max(lim, 0.05)
    im = ax.imshow(mat, aspect="auto", cmap="coolwarm", vmin=-lim, vmax=lim)
    ax.set_xticks(range(len(keep_res)))
    ax.set_xticklabels([str(r) for r in keep_res], fontsize=8)
    ax.set_yticks(range(len(muts)))
    ax.set_yticklabels(muts, fontsize=8)
    ax.set_xlabel("Protein residue (auth)")
    ax.set_title("Δ DOR–protein H-bond occupancy vs WT")
    fig.colorbar(im, ax=ax, fraction=0.03, label="Δ occupancy")
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def _plot_dccm_deltas(npy: Path, plots: Path) -> None:
    labels = [a for a in NNIBP_AUTH]
    for mut in STORY_MUTATIONS:
        path = npy / f"dccm_delta_{mut.replace('+', '_')}_minus_wt.npy"
        if not path.is_file():
            continue
        delta = np.load(path)
        fig, ax = plt.subplots(figsize=(5.2, 4.4))
        lim = float(np.max(np.abs(delta))) if delta.size else 1.0
        lim = max(lim, 0.05)
        im = ax.imshow(delta, cmap="coolwarm", vmin=-lim, vmax=lim)
        ax.set_xticks(range(len(labels)))
        ax.set_yticks(range(len(labels)))
        ax.set_xticklabels(labels, fontsize=7, rotation=90)
        ax.set_yticklabels(labels, fontsize=7)
        ax.set_title(f"NNIBP DCCM Δ  ({mut} − WT)")
        fig.colorbar(im, ax=ax, fraction=0.046)
        fig.tight_layout()
        fig.savefig(plots / f"dccm_delta_{mut.replace('+', '_')}.png", dpi=200)
        plt.close(fig)


def _plot_pose_pops(pop_g: pd.DataFrame, out: Path) -> None:
    muts = sorted(pop_g["mutation"].unique(), key=_mutation_sort_key)
    clusters = sorted(pop_g["cluster"].unique())
    fig, ax = plt.subplots(figsize=(11, 4.5))
    x = np.arange(len(muts))
    width = 0.8 / max(1, len(clusters))
    for i, c in enumerate(clusters):
        vals = []
        errs = []
        for mut in muts:
            row = pop_g[(pop_g["mutation"] == mut) & (pop_g["cluster"] == c)]
            if row.empty:
                vals.append(0.0)
                errs.append(0.0)
            else:
                vals.append(float(row["fraction_mean"].iloc[0]))
                errs.append(float(row["fraction_sem"].iloc[0]))
        ax.bar(x + i * width, vals, width=width, yerr=errs, capsize=1.5, label=f"C{c}")
    ax.set_xticks(x + width * (len(clusters) - 1) / 2)
    ax.set_xticklabels(muts, rotation=90, fontsize=8)
    ax.set_ylabel("Cluster fraction")
    ax.set_title("DOR pose-cluster populations (mean ± SEM)")
    ax.legend(fontsize=8, ncol=4)
    ax.grid(alpha=0.25, axis="y", linestyle=":")
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def _plot_network_deltas(net: pd.DataFrame, out: Path) -> None:
    focus = [m for m in STORY_MUTATIONS if m in set(net["mutation"])]
    if not focus:
        return
    # top |delta| edges per story mutant
    fig, axes = plt.subplots(1, len(focus), figsize=(3.6 * len(focus), 4.2), squeeze=False)
    for ax, mut in zip(axes[0], focus):
        sub = net[net["mutation"] == mut].copy()
        sub["abs"] = sub["delta_freq"].abs()
        top = sub.nlargest(8, "abs")
        labels = [f"{int(a)}–{int(b)}" for a, b in zip(top["auth_i"], top["auth_j"])]
        colors = ["#c0392b" if v > 0 else "#2c6fbb" for v in top["delta_freq"]]
        ax.barh(range(len(top)), top["delta_freq"], color=colors)
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(labels, fontsize=8)
        ax.axvline(0, color="0.4", lw=0.8)
        ax.set_title(mut, fontsize=10, fontweight="bold")
        ax.set_xlabel("Δ contact freq vs WT")
    fig.suptitle("NNIBP residue–residue contact remodeling", fontweight="bold")
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def _load_metas_from_tables(out: Path) -> list[dict]:
    inv = pd.read_csv(out / "tables" / "replicate_inventory.csv")
    rmsf = pd.read_csv(out / "tables" / "ligand_rmsf_per_atom.csv")
    hb_path = out / "tables" / "dor_hbond_occupancy_per_rep.csv"
    hb = pd.read_csv(hb_path) if hb_path.is_file() else pd.DataFrame()
    metas = []
    for _, row in inv.iterrows():
        m = row.to_dict()
        m["rmsf_rows"] = rmsf[(rmsf["mutation"] == row["mutation"]) & (rmsf["replicate"] == row["replicate"])].to_dict("records")
        if not hb.empty:
            m["hbond_rows"] = hb[(hb["mutation"] == row["mutation"]) & (hb["replicate"] == row["replicate"])].to_dict("records")
        else:
            m["hbond_rows"] = []
        metas.append(m)
    return metas


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    args = _parse_args()
    out = args.output_dir
    tables = out / "tables"
    plots = out / "plots"
    npy = out / "npy"
    cfg = out / "config"
    for d in (tables, plots, npy, cfg):
        d.mkdir(parents=True, exist_ok=True)

    if args.from_tables:
        metas = _load_metas_from_tables(out)
        _aggregate_and_plot(out, metas, n_pose_clusters=int(args.n_pose_clusters))
        LOGGER.info("Replotted from tables → %s", out)
        return 0

    run_df = collect_md_results(args.manifest)
    if args.mutations:
        keep = {_display_mutation(m) for m in args.mutations}
        run_df = run_df[run_df["mutation"].map(_display_mutation).isin(keep)].copy()
    jobs = _prepare_profile_jobs(run_df)
    # Require pbcfix
    filtered = []
    for job in jobs:
        traj = Path(job["trajectory"])
        if "pbcfix" not in traj.name:
            pbc = pbcfix_dcd_for(raw_analysis_dcd_for(traj) if traj.exists() else traj)
            if not pbc.exists():
                LOGGER.warning("skip %s rep%s: missing pbcfix", job["mutation"], job["replicate"])
                continue
            job["trajectory"] = str(pbc)
        filtered.append(job)
    jobs = filtered
    if not jobs:
        LOGGER.error("No jobs with pbcfix DCDs")
        return 1

    timing_rows = []
    for job in jobs:
        traj_path = Path(job["trajectory"])
        raw_dcd = raw_analysis_dcd_for(traj_path)
        timing_dcd = raw_dcd if raw_dcd.exists() else traj_path
        safe = str(job["safe_label"])
        rep = int(job["replicate"])
        call = infer_production_ns(
            dcd_path=timing_dcd,
            json_path=traj_path.parent / f"{safe}_rep{rep:02d}.json",
            state_csv_path=traj_path.parent / f"{safe}_rep{rep:02d}_md_state.csv",
            mutation=_display_mutation(job["mutation"]),
            replicate=rep,
        )
        job["production_ps"] = float(call.production_ns) * 1000.0
        job["timing_source"] = call.source
        timing_rows.append(call.__dict__)
    pd.DataFrame(timing_rows).to_csv(tables / "timing_audit.csv", index=False)

    metas: list[dict] = []
    failures: list[dict] = []
    LOGGER.info("Processing %d trajectories with %d workers", len(jobs), args.workers)
    with ProcessPoolExecutor(max_workers=max(1, int(args.workers))) as pool:
        futs = {
            pool.submit(
                _process_job,
                job,
                int(args.resid_offset),
                int(args.frame_stride),
                int(args.pocket_max_frames),
                int(args.hbond_max_frames),
                float(CONTACT_CUTOFF_A),
                str(npy),
            ): job
            for job in jobs
        }
        done = 0
        for fut in as_completed(futs):
            job = futs[fut]
            done += 1
            try:
                meta = fut.result()
                metas.append(meta)
                LOGGER.info(
                    "[%d/%d] %s rep%s  RMSF̄=%.2f Å  pocket=%.0f  Hbonds̄=%.1f",
                    done, len(jobs), meta["mutation"], meta["replicate"],
                    meta["rmsf_mean_angstrom"], meta["pocket_volume_mean"], meta["hbond_count_mean"],
                )
            except Exception as exc:
                failures.append(
                    {
                        "mutation": _display_mutation(job["mutation"]),
                        "replicate": int(job["replicate"]),
                        "error": str(exc),
                    }
                )
                LOGGER.warning("[%d/%d] FAILED %s rep%s: %s", done, len(jobs), job["mutation"], job["replicate"], exc)

    if failures:
        pd.DataFrame(failures).to_csv(tables / "failures.csv", index=False)
    if not metas:
        LOGGER.error("All jobs failed")
        return 1

    _aggregate_and_plot(out, metas, n_pose_clusters=int(args.n_pose_clusters))

    payload = {
        "manifest": str(args.manifest),
        "n_jobs": len(jobs),
        "n_ok": len(metas),
        "n_failed": len(failures),
        "resid_offset": args.resid_offset,
        "nnibp_auth": list(NNIBP_AUTH),
        "contact_cutoff_angstrom": CONTACT_CUTOFF_A,
        "n_pose_clusters": args.n_pose_clusters,
        "pocket_max_frames": args.pocket_max_frames,
        "hbond_max_frames": args.hbond_max_frames,
        "libraries": ["mdtraj", "MDAnalysis", "numpy", "scipy"],
        "notes": [
            "Uses *_analysis_pbcfix.dcd only.",
            "sklearn avoided (env pyarrow conflict); pose clusters via scipy kmeans2.",
            "PCA is SVD on per-replicate mean NNIBP Cα coordinates.",
            "Contact occupancy vs DOR is already in the draft; this suite adds H-bonds/RMSF/PCA/DCCM/networks.",
        ],
    }
    (cfg / "run_config.json").write_text(json.dumps(payload, indent=2) + "\n")
    LOGGER.info("Wrote %s (%d ok, %d failed)", out, len(metas), len(failures))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
