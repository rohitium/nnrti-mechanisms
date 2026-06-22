#!/usr/bin/env python3
"""Generate 100 ns triplet contact-story figures using pooled contact fractions and mean traces.

For each (WT, comparator, DRM) triplet this script:
1) Computes residue-level contact fractions across all frames and all replicates (pooled).
2) Picks a story residue where canonical contact is lower than WT/comparator.
3) Plots mean distance traces (average over all replicates) for the story residue.
4) Adds a raw contact heatmap (pooled contact fraction) for all contacted residues.
5) Labels legend with DOR fold-change values from DRM-susceptibilities sheet.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import textwrap
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from threading import Lock

import numpy as np
import pandas as pd

os.environ.setdefault("MPLBACKEND", "Agg")
_PLOT_LOCK = Lock()


@dataclass
class ReplicateMeta:
    mutation: str
    replicate: int
    output_json: Path
    topology_pdb: Path
    analysis_dcd: Path
    total_ns: float
    timing_source: str


def _steps_to_ns(steps: float | int | None, timestep_fs: float = 2.0) -> float:
    try:
        v = float(steps)
    except Exception:
        return np.nan
    if not np.isfinite(v):
        return np.nan
    return float(v * timestep_fs / 1_000_000.0)


def _resolve_local_path(path_like: str | Path | None, repo_root: Path) -> Path | None:
    if path_like is None:
        return None
    p = Path(str(path_like))
    if p.exists():
        return p
    marker = "nnrti-mechanisms/"
    text = str(p)
    if marker in text:
        mapped = repo_root / text.split(marker, 1)[1]
        if mapped.exists():
            return mapped
    rel = repo_root / str(p)
    if rel.exists():
        return rel
    return p


def _resolve_dcd_with_fallback(dcd_path: Path | None) -> Path | None:
    if dcd_path is None:
        return None
    if dcd_path.exists():
        return dcd_path
    name = dcd_path.name
    candidates: list[Path] = []
    if name.endswith("_analysis.dcd"):
        candidates.append(dcd_path.with_name(name.replace("_analysis.dcd", "_analysis.10ns.bak")))
    candidates.append(dcd_path.with_suffix(dcd_path.suffix + ".bak"))
    candidates.append(dcd_path.with_name(name + ".bak"))
    for c in candidates:
        if c.exists():
            return c
    return dcd_path


def _infer_total_ns_from_state_csv(path: Path | None) -> float:
    if path is None or not path.exists():
        return np.nan
    try:
        sdf = pd.read_csv(path)
    except Exception:
        return np.nan
    if sdf.empty:
        return np.nan
    step_col = None
    for c in ['#"Step"', '"#Step"', "Step"]:
        if c in sdf.columns:
            step_col = c
            break
    if step_col is None:
        return np.nan
    steps = pd.to_numeric(sdf[step_col], errors="coerce").dropna()
    if steps.empty:
        return np.nan
    return _steps_to_ns(float(steps.max()), timestep_fs=2.0)


def _load_replicate_meta(manifest_csv: Path, needed_mutations: set[str]) -> list[ReplicateMeta]:
    repo_root = Path(__file__).resolve().parents[3]
    mf = pd.read_csv(manifest_csv)
    req_cols = {"mutation", "replicate", "output_json"}
    missing = req_cols - set(mf.columns)
    if missing:
        raise ValueError(f"Manifest missing required columns: {sorted(missing)}")

    out: list[ReplicateMeta] = []
    for _, row in mf.sort_values(["mutation", "replicate"]).iterrows():
        mutation = str(row["mutation"])
        if mutation not in needed_mutations:
            continue
        replicate = int(pd.to_numeric(row["replicate"], errors="coerce"))
        out_json = _resolve_local_path(row["output_json"], repo_root=repo_root)
        if out_json is None or not out_json.exists():
            continue
        try:
            data = json.loads(out_json.read_text())
        except Exception:
            continue

        topo = _resolve_local_path(data.get("analysis_topology_pdb"), repo_root=repo_root)
        dcd = _resolve_local_path(data.get("analysis_dcd"), repo_root=repo_root)
        dcd = _resolve_dcd_with_fallback(dcd)
        if topo is None or dcd is None or (not topo.exists()) or (not dcd.exists()):
            continue

        ns_json = _steps_to_ns(data.get("md_production_steps_completed", data.get("md_production_steps")), timestep_fs=2.0)
        state_csv = _resolve_local_path(data.get("state_csv"), repo_root=repo_root)
        ns_state = _infer_total_ns_from_state_csv(state_csv)

        has_state = bool(np.isfinite(ns_state) and ns_state > 0)
        has_json = bool(np.isfinite(ns_json) and ns_json > 0)
        if has_state and has_json:
            # Prefer the longer duration when sources disagree (for example,
            # stale state CSV after resumed production that completed in JSON).
            if float(ns_state) >= float(ns_json):
                total_ns = float(ns_state)
                timing_source = "state_csv"
            else:
                total_ns = float(ns_json)
                timing_source = "json_steps_gt_state_csv"
        elif has_state:
            total_ns = float(ns_state)
            timing_source = "state_csv"
        elif has_json:
            total_ns = float(ns_json)
            timing_source = "json_steps"
        else:
            total_ns = np.nan
            timing_source = "unknown"

        out.append(
            ReplicateMeta(
                mutation=mutation,
                replicate=replicate,
                output_json=out_json,
                topology_pdb=topo,
                analysis_dcd=dcd,
                total_ns=total_ns,
                timing_source=timing_source,
            )
        )
    return out


def _parse_triplets(text: str) -> list[tuple[str, str, str]]:
    triplets: list[tuple[str, str, str]] = []
    for block in str(text).split(";"):
        block = block.strip()
        if not block:
            continue
        toks = [x.strip() for x in block.split(",") if x.strip()]
        if len(toks) != 3:
            raise ValueError(f"Triplet must have 3 comma-separated mutations: {block}")
        triplets.append((toks[0], toks[1], toks[2]))
    if not triplets:
        raise ValueError("No triplets parsed.")
    return triplets


def _parse_int_csv(text: str | None) -> list[int]:
    if text is None:
        return []
    out: list[int] = []
    for tok in re.split(r"[,\s;]+", str(text).strip()):
        if not tok:
            continue
        out.append(int(tok))
    return out


def _normalize_mutation_token(text: str) -> str:
    t = str(text).strip().upper()
    if not t or t == "NAN":
        return ""
    t = t.replace(" ", "")
    t = t.replace(",", "+")
    t = re.sub(r"\++", "+", t)
    return t


def _load_dor_fold_map(xlsx_path: Path) -> dict[str, float]:
    # Sheet has a blank first row; actual headers start on row 2.
    df = pd.read_excel(xlsx_path, header=1)
    if df.shape[1] < 3:
        raise ValueError(f"Unexpected susceptibility sheet format: {xlsx_path}")
    df = df.iloc[:, :3].copy()
    df.columns = ["mutation_raw", "rpv_fold", "dor_fold"]
    df["mutation"] = df["mutation_raw"].map(_normalize_mutation_token)
    df["dor_fold"] = pd.to_numeric(df["dor_fold"], errors="coerce")
    df = df[(df["mutation"] != "") & df["dor_fold"].notna()].copy()
    out = {str(m): float(v) for m, v in zip(df["mutation"], df["dor_fold"]) if np.isfinite(v)}
    out["WT"] = 1.0
    return out


def _residue_label(auth_resid: int, resname: str | None) -> str:
    if resname and str(resname).strip():
        return f"{str(resname).strip()}{int(auth_resid)}"
    return f"Res{int(auth_resid)}"


def _position_label(auth_resid: int) -> str:
    return f"Pos{int(auth_resid)}"


def _aa1_to_aa3(code: str) -> str:
    m = {
        "A": "ALA",
        "R": "ARG",
        "N": "ASN",
        "D": "ASP",
        "C": "CYS",
        "Q": "GLN",
        "E": "GLU",
        "G": "GLY",
        "H": "HIS",
        "I": "ILE",
        "L": "LEU",
        "K": "LYS",
        "M": "MET",
        "F": "PHE",
        "P": "PRO",
        "S": "SER",
        "T": "THR",
        "W": "TRP",
        "Y": "TYR",
        "V": "VAL",
    }
    return m.get(str(code).upper(), "UNK")


def _parse_mutation_ops(mutation: str) -> list[tuple[str, int, str]]:
    ops: list[tuple[str, int, str]] = []
    token = _normalize_mutation_token(mutation)
    if token in {"", "WT"}:
        return ops
    for part in token.split("+"):
        mt = re.match(r"^([A-Z])(\d+)([A-Z])$", str(part).strip().upper())
        if mt is None:
            continue
        ops.append((str(mt.group(1)), int(mt.group(2)), str(mt.group(3))))
    return ops


def _infer_resname_map_for_position(
    triplet: tuple[str, str, str],
    auth_resid: int,
) -> tuple[dict[str, str], bool]:
    pos = int(auth_resid)
    ref_aa1 = ""
    for mut in triplet:
        for src, mpos, _dst in _parse_mutation_ops(mut):
            if int(mpos) == pos:
                ref_aa1 = str(src).upper()
                break
        if ref_aa1:
            break

    out: dict[str, str] = {}
    has_pos_mutation = False
    for mut in triplet:
        mut_ops = _parse_mutation_ops(mut)
        dest = ""
        for _src, mpos, dst in mut_ops:
            if int(mpos) == pos:
                dest = str(dst).upper()
                has_pos_mutation = True
                break
        if dest:
            out[str(mut)] = _aa1_to_aa3(dest)
        elif ref_aa1:
            out[str(mut)] = _aa1_to_aa3(ref_aa1)
        else:
            out[str(mut)] = "UNK"
    return out, bool(has_pos_mutation)


def _label_resname_for_position(
    mut_occ: pd.DataFrame,
    wt_mutation: str,
    traj_resid: int,
    auth_resid: int,
) -> str:
    sel_wt = mut_occ[
        (mut_occ["mutation"] == str(wt_mutation))
        & (pd.to_numeric(mut_occ["traj_resid"], errors="coerce") == int(traj_resid))
        & (pd.to_numeric(mut_occ["auth_resid"], errors="coerce") == int(auth_resid))
    ].copy()
    source = sel_wt
    if source.empty:
        source = mut_occ[
            (pd.to_numeric(mut_occ["traj_resid"], errors="coerce") == int(traj_resid))
            & (pd.to_numeric(mut_occ["auth_resid"], errors="coerce") == int(auth_resid))
        ].copy()
    if source.empty:
        return ""
    names = source["resname"].astype(str).str.strip()
    names = names[(names != "") & (names.str.upper() != "NAN")]
    if names.empty:
        return ""
    vc = names.value_counts(dropna=False)
    return str(vc.index[0]).upper()


def _display_resname_for_position(
    mut_occ: pd.DataFrame,
    triplet: tuple[str, str, str],
    wt_mutation: str,
    traj_resid: int,
    auth_resid: int,
) -> str:
    res_map, has_pos_mut = _infer_resname_map_for_position(triplet=triplet, auth_resid=int(auth_resid))
    if has_pos_mut:
        wt_name = str(res_map.get(str(wt_mutation), "")).strip().upper()
        if wt_name and wt_name != "UNK":
            return wt_name
    return _label_resname_for_position(
        mut_occ=mut_occ,
        wt_mutation=str(wt_mutation),
        traj_resid=int(traj_resid),
        auth_resid=int(auth_resid),
    )


def _select_residue_sidechain(universe, traj_resid: int):
    residues = universe.select_atoms(f"protein and resid {int(traj_resid)}").residues
    if len(residues) == 0:
        return universe.select_atoms(f"protein and resid {int(traj_resid)} and not name H*")

    picked: list[int] = []
    for res in residues:
        rag = res.atoms
        hv = rag.select_atoms("not name H*")
        picked.extend([int(i) for i in hv.indices.tolist()])

    if not picked:
        return universe.select_atoms(f"protein and resid {int(traj_resid)} and not name H*")
    uniq = np.asarray(sorted(set(picked)), dtype=int)
    return universe.atoms[uniq]


def _compute_triplet_contact_stats(
    metas: list[ReplicateMeta],
    mutation_triplet: tuple[str, str, str],
    ligand_resname: str,
    resid_offset: int,
    contact_cutoff: float,
    window_ns: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    import MDAnalysis as mda
    from MDAnalysis import transformations as trans
    from MDAnalysis.lib.distances import capped_distance

    muts = set(mutation_triplet)
    meta_sel = [m for m in metas if m.mutation in muts]
    if not meta_sel:
        raise ValueError(f"No valid replicates found for triplet: {mutation_triplet}")

    rep_counts: dict[tuple[str, int, int], int] = {}
    rep_frames: dict[tuple[str, int], int] = {}
    res_meta: dict[int, str] = {}
    timing_rows: list[dict[str, object]] = []

    for m in sorted(meta_sel, key=lambda x: (x.mutation, x.replicate)):
        u = mda.Universe(str(m.topology_pdb), str(m.analysis_dcd), format="DCD")
        lig = u.select_atoms(f"resname {ligand_resname} and not name H*")
        prot = u.select_atoms("protein and not name H*")
        if lig.n_atoms == 0 or prot.n_atoms == 0:
            continue

        try:
            anchor = u.select_atoms("protein")
            if anchor.n_atoms == 0:
                anchor = u.atoms
            u.trajectory.add_transformations(
                trans.NoJump(check_continuity=False),
                trans.center_in_box(anchor, center="geometry", wrap=False),
            )
        except Exception:
            pass

        n_frames = len(u.trajectory)
        if n_frames < 2:
            continue
        total_ns = float(m.total_ns) if np.isfinite(m.total_ns) and m.total_ns > 0 else float(window_ns)
        t_ns = np.linspace(0.0, total_ns, n_frames)
        keep_idx = np.where(t_ns <= float(window_ns))[0]
        if keep_idx.size < 2:
            keep_idx = np.arange(n_frames, dtype=int)

        atom_to_resid = np.asarray([int(a.resid) for a in prot.atoms], dtype=int)
        rep_key = (m.mutation, int(m.replicate))
        rep_frames[rep_key] = int(len(keep_idx))
        timing_rows.append(
            {
                "mutation": m.mutation,
                "replicate": int(m.replicate),
                "n_frames_total": int(n_frames),
                "n_frames_window": int(len(keep_idx)),
                "total_ns_used": float(total_ns),
                "timing_source": m.timing_source,
                "analysis_dcd": str(m.analysis_dcd),
            }
        )

        for fi in keep_idx.tolist():
            u.trajectory[int(fi)]
            pairs = capped_distance(
                prot.positions,
                lig.positions,
                max_cutoff=float(contact_cutoff),
                min_cutoff=0.0,
                box=u.dimensions,
                return_distances=False,
            )
            if pairs is None or len(pairs) == 0:
                continue
            touched = set(atom_to_resid[np.asarray(pairs)[:, 0]].tolist())
            for tr in touched:
                k = (m.mutation, int(m.replicate), int(tr))
                rep_counts[k] = rep_counts.get(k, 0) + 1
                if int(tr) not in res_meta:
                    ag = u.select_atoms(f"protein and resid {int(tr)}")
                    res_meta[int(tr)] = str(ag.residues[0].resname) if ag.n_atoms > 0 and len(ag.residues) else ""

    if not rep_frames:
        raise ValueError(f"No usable frames for triplet {mutation_triplet}")

    rep_rows: list[dict[str, object]] = []
    for (mutation, replicate), nfr in sorted(rep_frames.items()):
        touched_resids = sorted({k[2] for k in rep_counts.keys() if k[0] == mutation and k[1] == replicate})
        for tr in touched_resids:
            cnt = int(rep_counts.get((mutation, replicate, tr), 0))
            auth = int(tr) - int(resid_offset)
            rep_rows.append(
                {
                    "mutation": mutation,
                    "replicate": int(replicate),
                    "traj_resid": int(tr),
                    "auth_resid": int(auth),
                    "resname": str(res_meta.get(int(tr), "")),
                    "n_contact_frames": cnt,
                    "n_total_frames": int(nfr),
                    "occupancy": float(cnt / max(1, nfr)),
                }
            )
    rep_occ = pd.DataFrame(rep_rows)
    if rep_occ.empty:
        raise ValueError(f"No residue contacts detected for triplet {mutation_triplet}")

    # Both unweighted replicate mean and pooled occupancy (all frames across reps).
    mut_occ = (
        rep_occ.groupby(["mutation", "traj_resid", "auth_resid", "resname"], as_index=False)
        .agg(
            occupancy_mean=("occupancy", "mean"),
            occupancy_std=("occupancy", "std"),
            n_replicates=("replicate", "nunique"),
            n_contact_frames_total=("n_contact_frames", "sum"),
            n_total_frames_total=("n_total_frames", "sum"),
        )
    )
    mut_occ["occupancy_sem"] = mut_occ["occupancy_std"] / np.sqrt(mut_occ["n_replicates"].clip(lower=1))
    mut_occ["occupancy_pooled"] = (
        pd.to_numeric(mut_occ["n_contact_frames_total"], errors="coerce")
        / pd.to_numeric(mut_occ["n_total_frames_total"], errors="coerce").clip(lower=1)
    )
    timing_df = pd.DataFrame(timing_rows)
    return rep_occ, mut_occ, timing_df


def _choose_story_residue(
    mut_occ: pd.DataFrame,
    triplet: tuple[str, str, str],
    min_wt_pooled: float,
    min_neg_pooled: float,
    max_can_pooled: float,
    min_any_mean_occ: float,
) -> pd.Series:
    wt, neg, can = triplet
    piv = (
        mut_occ.pivot_table(
            index=["traj_resid", "auth_resid"],
            columns="mutation",
            values="occupancy_pooled",
            aggfunc="mean",
        )
        .reset_index()
    )
    piv_mean = (
        mut_occ.pivot_table(
            index=["traj_resid", "auth_resid"],
            columns="mutation",
            values="occupancy_mean",
            aggfunc="mean",
        )
        .reset_index()
    )
    for c in [wt, neg, can]:
        if c not in piv.columns:
            piv[c] = np.nan
        if c not in piv_mean.columns:
            piv_mean[c] = np.nan
    piv_mean["max_any_mean_occ"] = piv_mean[[wt, neg, can]].max(axis=1, skipna=True)
    piv = piv.merge(
        piv_mean[["traj_resid", "auth_resid", "max_any_mean_occ"]],
        on=["traj_resid", "auth_resid"],
        how="left",
    )
    wt_vals = pd.to_numeric(piv[wt], errors="coerce")
    neg_vals = pd.to_numeric(piv[neg], errors="coerce")
    can_vals = pd.to_numeric(piv[can], errors="coerce")

    # Directional story scores:
    #   canonical_lower  = canonical contact loss vs both WT and comparator
    #   canonical_higher = canonical contact gain vs both WT and comparator
    piv["down_score"] = np.minimum(wt_vals, neg_vals) - can_vals
    piv["up_score"] = can_vals - np.maximum(wt_vals, neg_vals)
    piv["score"] = np.maximum(pd.to_numeric(piv["down_score"], errors="coerce"), pd.to_numeric(piv["up_score"], errors="coerce"))
    piv["direction"] = np.where(
        pd.to_numeric(piv["down_score"], errors="coerce") >= pd.to_numeric(piv["up_score"], errors="coerce"),
        "canonical_lower",
        "canonical_higher",
    )
    piv = piv.replace([np.inf, -np.inf], np.nan).dropna(
        subset=[wt, neg, can, "down_score", "up_score", "score", "max_any_mean_occ"]
    ).copy()
    if piv.empty:
        raise ValueError(f"No comparable residues across triplet {triplet}")

    wt_vals = pd.to_numeric(piv[wt], errors="coerce")
    neg_vals = pd.to_numeric(piv[neg], errors="coerce")
    can_vals = pd.to_numeric(piv[can], errors="coerce")

    down_ok = (
        (wt_vals >= float(min_wt_pooled))
        & (neg_vals >= float(min_neg_pooled))
        & (can_vals <= float(max_can_pooled))
        & (pd.to_numeric(piv["down_score"], errors="coerce") > 0.0)
    )
    up_ok = (
        (can_vals >= float(min_wt_pooled))
        & (wt_vals <= float(max_can_pooled))
        & (neg_vals <= float(max_can_pooled))
        & (pd.to_numeric(piv["up_score"], errors="coerce") > 0.0)
    )
    constrained = piv[
        (pd.to_numeric(piv["max_any_mean_occ"], errors="coerce") > float(min_any_mean_occ))
        & (down_ok | up_ok)
    ].copy()
    if not constrained.empty:
        return constrained.sort_values(["score", "max_any_mean_occ", wt, neg], ascending=[False, False, False, False]).iloc[0]

    fallback = piv[pd.to_numeric(piv["max_any_mean_occ"], errors="coerce") > float(min_any_mean_occ)].copy()
    if not fallback.empty:
        return fallback.sort_values(["score", "max_any_mean_occ", wt, neg], ascending=[False, False, False, False]).iloc[0]
    return piv.sort_values(["score", wt, neg], ascending=[False, False, False]).iloc[0]


def _extract_distance_trace(
    meta: ReplicateMeta,
    auth_resid: int,
    ligand_resname: str,
    resid_offset: int,
    window_ns: float,
) -> tuple[np.ndarray, np.ndarray]:
    import MDAnalysis as mda
    from MDAnalysis import transformations as trans
    from MDAnalysis.lib.distances import distance_array

    u = mda.Universe(str(meta.topology_pdb), str(meta.analysis_dcd), format="DCD")
    lig = u.select_atoms(f"resname {ligand_resname} and not name H*")
    if lig.n_atoms == 0:
        raise ValueError(f"Empty ligand selection for {meta.mutation} rep{meta.replicate}")

    traj_resid = int(auth_resid) + int(resid_offset)
    sc = _select_residue_sidechain(u, traj_resid=traj_resid)
    if sc.n_atoms == 0:
        raise ValueError(f"No residue atoms for auth {auth_resid} (traj {traj_resid}) in {meta.mutation}")

    try:
        anchor = u.select_atoms("protein")
        if anchor.n_atoms == 0:
            anchor = u.atoms
        u.trajectory.add_transformations(
            trans.NoJump(check_continuity=False),
            trans.center_in_box(anchor, center="geometry", wrap=False),
        )
    except Exception:
        pass

    n_frames = len(u.trajectory)
    if n_frames < 2:
        raise ValueError(f"Too few frames in {meta.analysis_dcd}")

    total_ns = float(meta.total_ns) if np.isfinite(meta.total_ns) and meta.total_ns > 0 else float(window_ns)
    t_ns = np.linspace(0.0, total_ns, n_frames)
    keep = t_ns <= float(window_ns)
    if int(np.sum(keep)) < 2:
        keep = np.ones(n_frames, dtype=bool)

    dvals = np.full(n_frames, np.nan, dtype=float)
    for i, _ in enumerate(u.trajectory):
        dvals[i] = float(distance_array(sc.positions, lig.positions, box=u.dimensions).min())
    return t_ns[keep], dvals[keep]


def _interp_trace_to_grid(x: np.ndarray, y: np.ndarray, grid: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    x = x[keep]
    y = y[keep]
    if x.size < 2:
        return np.full_like(grid, np.nan, dtype=float)
    mono = np.r_[True, np.diff(x) > 0]
    x = x[mono]
    y = y[mono]
    if x.size < 2:
        return np.full_like(grid, np.nan, dtype=float)
    yi = np.interp(grid, x, y, left=np.nan, right=np.nan)
    yi[(grid < x.min()) | (grid > x.max())] = np.nan
    return yi


def _extract_mutation_mean_trace(
    metas: list[ReplicateMeta],
    mutation: str,
    auth_resid: int,
    ligand_resname: str,
    resid_offset: int,
    window_ns: float,
    n_grid: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[dict[str, object]]]:
    reps = [m for m in metas if m.mutation == mutation]
    if not reps:
        raise ValueError(f"No replicate metadata for mutation: {mutation}")

    grid = np.linspace(0.0, float(window_ns), int(n_grid))
    ys: list[np.ndarray] = []
    rep_rows: list[dict[str, object]] = []
    for m in sorted(reps, key=lambda x: x.replicate):
        x, dist_by_res = _extract_rep_distance_traces_mdtraj(
            meta=m,
            auth_resids=[int(auth_resid)],
            ligand_resname=str(ligand_resname),
            resid_offset=int(resid_offset),
            window_ns=float(window_ns),
        )
        y = np.asarray(dist_by_res.get(int(auth_resid), np.full_like(x, np.nan, dtype=float)), dtype=float)
        yi = _interp_trace_to_grid(x, y, grid)
        ys.append(yi)
        rep_rows.append(
            {
                "mutation": mutation,
                "replicate": int(m.replicate),
                "n_points": int(np.isfinite(yi).sum()),
                "timing_source": m.timing_source,
            }
        )
    stack = np.vstack(ys)
    n_eff = np.isfinite(stack).sum(axis=0)
    sum_vals = np.nansum(stack, axis=0)
    mean = np.divide(sum_vals, n_eff, out=np.full_like(sum_vals, np.nan, dtype=float), where=n_eff > 0)

    centered = stack - mean[None, :]
    centered[~np.isfinite(stack)] = np.nan
    sq = np.nansum(centered ** 2, axis=0)
    sd = np.divide(sq, np.maximum(n_eff - 1, 1), out=np.full_like(sq, np.nan, dtype=float), where=n_eff > 1) ** 0.5
    sem = np.divide(sd, np.sqrt(n_eff), out=np.full_like(sd, np.nan), where=n_eff > 1)
    return grid, mean, sem, rep_rows


def _is_hydrogen_atom_mdtraj(atom) -> bool:
    if getattr(atom, "element", None) is not None:
        symbol = str(getattr(atom.element, "symbol", "")).strip().upper()
        if symbol == "H":
            return True
    return str(getattr(atom, "name", "")).strip().upper().startswith("H")


def _select_residue_sidechain_indices_mdtraj(topology, traj_resid: int) -> np.ndarray:
    traj_resid = int(traj_resid)
    residues = [r for r in topology.residues if r.is_protein and int(r.resSeq) == traj_resid]
    if not residues:
        return np.asarray([], dtype=int)

    picked: list[int] = []
    for res in residues:
        atoms = list(res.atoms)
        heavy = [a.index for a in atoms if not _is_hydrogen_atom_mdtraj(a)]
        picked.extend(heavy)
    return np.asarray(sorted(set(int(i) for i in picked)), dtype=int)


def _extract_rep_distance_traces_mdtraj(
    meta: ReplicateMeta,
    auth_resids: list[int],
    ligand_resname: str,
    resid_offset: int,
    window_ns: float,
) -> tuple[np.ndarray, dict[int, np.ndarray]]:
    import mdtraj as md

    traj = md.load_dcd(str(meta.analysis_dcd), top=str(meta.topology_pdb))
    n_frames = int(traj.n_frames)
    if n_frames < 2:
        raise ValueError(f"Too few frames in {meta.analysis_dcd}")

    total_ns = float(meta.total_ns) if np.isfinite(meta.total_ns) and meta.total_ns > 0 else float(window_ns)
    t_ns = np.linspace(0.0, total_ns, n_frames)
    keep = t_ns <= float(window_ns)
    if int(np.sum(keep)) < 2:
        keep = np.ones(n_frames, dtype=bool)
    traj = traj[keep]
    t_sel = t_ns[keep]

    try:
        # Reconstruct molecule connectivity directly from bonds (robust against
        # image_molecules splitting protein/ligand into different images).
        traj.make_molecules_whole(inplace=True)
    except Exception:
        pass

    top = traj.topology
    lig_idx = np.asarray(
        [
            a.index
            for a in top.atoms
            if str(a.residue.name).strip() == str(ligand_resname) and not _is_hydrogen_atom_mdtraj(a)
        ],
        dtype=int,
    )
    if lig_idx.size == 0:
        raise ValueError(f"Empty ligand selection ({ligand_resname}) for {meta.mutation} rep{meta.replicate}")

    prot_idx = np.asarray(
        [a.index for a in top.atoms if bool(getattr(a.residue, "is_protein", False))],
        dtype=int,
    )
    if (
        lig_idx.size > 0
        and prot_idx.size > 0
        and getattr(traj, "unitcell_lengths", None) is not None
    ):
        # Per-frame nearest-image placement of ligand around protein COM.
        for fi in range(traj.n_frames):
            box = traj.unitcell_lengths[fi]
            if box is None or not np.all(np.isfinite(box)) or not np.all(box > 0):
                continue
            prot_c = traj.xyz[fi, prot_idx].mean(axis=0)
            lig_c = traj.xyz[fi, lig_idx].mean(axis=0)
            delta = lig_c - prot_c
            traj.xyz[fi, lig_idx] += -box * np.round(delta / box)

    out: dict[int, np.ndarray] = {}
    for auth_resid in auth_resids:
        traj_resid = int(auth_resid) + int(resid_offset)
        sc_idx = _select_residue_sidechain_indices_mdtraj(top, traj_resid=traj_resid)
        if sc_idx.size == 0:
            out[int(auth_resid)] = np.full(traj.n_frames, np.nan, dtype=float)
            continue
        pairs = np.asarray([(int(i), int(j)) for i in sc_idx for j in lig_idx], dtype=int)
        if pairs.size == 0:
            out[int(auth_resid)] = np.full(traj.n_frames, np.nan, dtype=float)
            continue
        # Use minimum-image distances; this is robust to residual image hops
        # that can survive trajectory reimaging in some replicates.
        d_nm = md.compute_distances(traj, pairs, periodic=True, opt=True)
        out[int(auth_resid)] = np.nanmin(d_nm, axis=1) * 10.0
    return t_sel, out


def _extract_mutation_mean_traces_mdtraj(
    metas: list[ReplicateMeta],
    mutation: str,
    auth_resids: list[int],
    ligand_resname: str,
    resid_offset: int,
    window_ns: float,
    n_grid: int,
) -> dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    reps = [m for m in metas if m.mutation == mutation]
    if not reps:
        raise ValueError(f"No replicate metadata for mutation: {mutation}")
    if not auth_resids:
        raise ValueError("No auth_resids provided for trace extraction.")

    grid = np.linspace(0.0, float(window_ns), int(n_grid))
    per_residue: dict[int, list[np.ndarray]] = {int(r): [] for r in auth_resids}

    for m in sorted(reps, key=lambda x: x.replicate):
        t_sel, d_by_res = _extract_rep_distance_traces_mdtraj(
            meta=m,
            auth_resids=[int(r) for r in auth_resids],
            ligand_resname=ligand_resname,
            resid_offset=int(resid_offset),
            window_ns=float(window_ns),
        )
        for auth_resid in auth_resids:
            y = np.asarray(d_by_res.get(int(auth_resid), np.full_like(t_sel, np.nan, dtype=float)), dtype=float)
            yi = _interp_trace_to_grid(t_sel, y, grid)
            per_residue[int(auth_resid)].append(yi)

    out: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for auth_resid in auth_resids:
        stack = np.vstack(per_residue[int(auth_resid)])
        n_eff = np.isfinite(stack).sum(axis=0)
        sum_vals = np.nansum(stack, axis=0)
        mean = np.divide(sum_vals, n_eff, out=np.full_like(sum_vals, np.nan, dtype=float), where=n_eff > 0)

        centered = stack - mean[None, :]
        centered[~np.isfinite(stack)] = np.nan
        sq = np.nansum(centered ** 2, axis=0)
        sd = np.divide(sq, np.maximum(n_eff - 1, 1), out=np.full_like(sq, np.nan, dtype=float), where=n_eff > 1) ** 0.5
        sem = np.divide(sd, np.sqrt(n_eff), out=np.full_like(sd, np.nan), where=n_eff > 1)
        out[int(auth_resid)] = (grid, mean, sem)
    return out


def _format_fold_label(mutation: str, fold_map: dict[str, float]) -> str:
    key = _normalize_mutation_token(mutation)
    v = fold_map.get(key, np.nan)
    if np.isfinite(v):
        if abs(v - round(v)) < 1e-8:
            return f"{mutation} ({int(round(v))}x)"
        return f"{mutation} ({v:.1f}x)"
    return mutation


def _wt_contact_region(auth_resid: int) -> tuple[int, str, str]:
    auth_resid = int(auth_resid)
    region_specs = [
        (0, "β6 strand", "#4c78a8", {95, 97}),
        (1, "Pocket entrance", "#f58518", {100, 101, 103, 179, 181}),
        (2, "103-108 loop", "#eeca3b", {102, 104, 105, 106, 107, 108}),
        (3, "Hydrophobic tunnel", "#54a24b", {188, 227, 229, 234}),
        (4, "β9-β10 hairpin", "#e45756", {180, 189, 190}),
        (5, "β12-β13 primer grip", "#b279a2", {223, 225, 228, 235, 236, 237}),
        (6, "Distal pocket wall", "#72b7b2", {318}),
    ]
    for order, label, color, members in region_specs:
        if auth_resid in members:
            return int(order), str(label), str(color)
    return 99, "other pocket contacts", "#9d9da1"


def _plot_triplet_figure(
    triplet: tuple[str, str, str],
    mut_occ: pd.DataFrame,
    story_row: pd.Series,
    mean_traces: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]],
    fold_map: dict[str, float],
    contact_cutoff: float,
    min_any_mean_occ_display: float,
    output_png: Path,
) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm

    wt, neg, can = triplet
    mut_order = [wt, neg, can]
    colors = {wt: "#4d4d4d", neg: "#1f77b4", can: "#d62728"}

    disp_filter = (
        mut_occ[mut_occ["mutation"].isin(mut_order)]
        .groupby(["traj_resid", "auth_resid"], as_index=False)["occupancy_mean"]
        .max()
    )
    disp_filter["occupancy_mean"] = pd.to_numeric(disp_filter["occupancy_mean"], errors="coerce")
    disp_filter = disp_filter[disp_filter["occupancy_mean"] > float(min_any_mean_occ_display)].copy()
    if disp_filter.empty:
        raise ValueError(
            f"No residues pass contact_mean > {min_any_mean_occ_display:.2f} for heatmap display in triplet {triplet}"
        )

    all_res = (
        disp_filter[["traj_resid", "auth_resid"]]
        .drop_duplicates()
        .sort_values(["auth_resid", "traj_resid"])
        .copy()
    )
    all_res["resname"] = all_res.apply(
        lambda r: _display_resname_for_position(
            mut_occ=mut_occ,
            triplet=triplet,
            wt_mutation=wt,
            traj_resid=int(r["traj_resid"]),
            auth_resid=int(r["auth_resid"]),
        ),
        axis=1,
    )
    all_res["label"] = all_res.apply(
        lambda r: _residue_label(int(r["auth_resid"]), str(r["resname"])),
        axis=1,
    )
    col_order = all_res["label"].tolist()
    key_lookup = all_res.copy()
    key_lookup["key"] = list(zip(key_lookup["traj_resid"], key_lookup["auth_resid"]))
    map_lab = dict(zip(key_lookup["key"], key_lookup["label"]))
    tmp = mut_occ.copy()
    tmp["key"] = list(zip(tmp["traj_resid"], tmp["auth_resid"]))
    tmp = tmp[tmp["key"].isin(set(map_lab.keys()))].copy()
    tmp["label"] = tmp["key"].map(map_lab)
    hm = (
        tmp.pivot_table(index="mutation", columns="label", values="occupancy_pooled", aggfunc="mean")
        .reindex(index=mut_order, columns=col_order)
    )

    fig, (ax_top, ax_bot) = plt.subplots(
        2,
        1,
        figsize=(13.8, 8.8),
        gridspec_kw={"height_ratios": [1.35, 1.0]},
        constrained_layout=True,
    )

    xmax = 0.0
    for mut in mut_order:
        x, y, sem = mean_traces[mut]
        xmax = max(xmax, float(np.nanmax(x)) if len(x) else 0.0)
        ax_top.plot(
            x,
            y,
            color=colors[mut],
            linewidth=2.1,
            alpha=0.95,
            label=_format_fold_label(mut, fold_map=fold_map),
        )
        lo = y - sem
        hi = y + sem
        ok = np.isfinite(lo) & np.isfinite(hi)
        if np.any(ok):
            ax_top.fill_between(x[ok], lo[ok], hi[ok], color=colors[mut], alpha=0.16, linewidth=0)

    ax_top.axhline(float(contact_cutoff), color="#666666", linestyle=":", linewidth=1.1, label=f"{contact_cutoff:.1f} A contact cutoff")
    ax_top.set_xlim(0.0, xmax if xmax > 0 else 100.0)
    ax_top.set_xlabel("Time (ns)")
    ax_top.set_ylabel("Min residue-DOR distance (A)")
    story_auth = int(story_row["auth_resid"])
    story_traj = int(story_row["traj_resid"])
    story_resname = _display_resname_for_position(
        mut_occ=mut_occ,
        triplet=triplet,
        wt_mutation=wt,
        traj_resid=story_traj,
        auth_resid=story_auth,
    )
    story_label = _residue_label(story_auth, story_resname)
    ax_top.set_title(story_label)
    ax_top.grid(alpha=0.22, linestyle=":")
    ax_top.legend(loc="upper right", frameon=True, fontsize=8)

    arr = hm.to_numpy(dtype=float)
    cmap = plt.get_cmap("cividis")
    boundaries = np.arange(0.5, 1.0001, 0.025)
    norm = BoundaryNorm(boundaries=boundaries, ncolors=cmap.N, clip=True)
    im = ax_bot.imshow(arr, aspect="auto", cmap=cmap, norm=norm)
    ax_bot.set_yticks(np.arange(len(mut_order)), [_format_fold_label(m, fold_map=fold_map) for m in mut_order])
    ax_bot.set_xticks(np.arange(len(col_order)), col_order, rotation=45, ha="right")
    ax_bot.set_xlabel("Residue")
    ax_bot.set_ylabel("Mutation")
    ax_bot.set_title("Mean Occupancy Heatmap")
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            v = arr[i, j]
            if not np.isfinite(v):
                continue
            txt_color = "white" if float(v) < 0.72 else "black"
            ax_bot.text(j, i, f"{float(v):.2f}", ha="center", va="center", fontsize=7, color=txt_color)
    cbar = fig.colorbar(im, ax=ax_bot, fraction=0.046, pad=0.02)
    cbar.set_ticks(np.arange(0.5, 1.01, 0.05))

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=320, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output_png}")


def _plot_wt_contacts_figure(
    wt_mutation: str,
    mut_occ: pd.DataFrame,
    metas: list[ReplicateMeta],
    ligand_resname: str,
    resid_offset: int,
    window_ns: float,
    n_grid: int,
    contact_cutoff: float,
    min_wt_mean_occ: float,
    wt_trace_auth_resids: list[int] | None,
    output_png: Path,
) -> None:
    import matplotlib.pyplot as plt

    wt = str(wt_mutation)
    wt_occ_all = mut_occ[mut_occ["mutation"] == wt].copy()
    wt_occ_all["occupancy_mean"] = pd.to_numeric(wt_occ_all["occupancy_mean"], errors="coerce")
    wt_occ_all["occupancy_sem"] = pd.to_numeric(wt_occ_all["occupancy_sem"], errors="coerce")
    wt_occ_all = wt_occ_all[np.isfinite(wt_occ_all["occupancy_mean"])].copy()
    if wt_occ_all.empty:
        raise ValueError(f"No contacted residues available for mutation {wt}")

    wt_occ_all[["region_order", "pocket_region", "region_color"]] = wt_occ_all["auth_resid"].apply(
        lambda v: pd.Series(_wt_contact_region(int(v)))
    )
    wt_occ_all = wt_occ_all.sort_values(["region_order", "auth_resid", "occupancy_mean"], ascending=[True, True, False]).reset_index(drop=True)
    wt_occ_all["label"] = wt_occ_all.apply(
        lambda r: _residue_label(int(r["auth_resid"]), str(r.get("resname", ""))),
        axis=1,
    )

    wt_occ_trace = wt_occ_all[wt_occ_all["occupancy_mean"] >= float(min_wt_mean_occ)].copy()
    if wt_occ_trace.empty:
        raise ValueError(f"No residues satisfy contact_mean >= {min_wt_mean_occ:.2f} for mutation {wt}")

    allowed_auth_all = {int(v) for v in wt_occ_all["auth_resid"].tolist()}
    requested_auth = [int(v) for v in (wt_trace_auth_resids or [])]
    if requested_auth:
        auth_resids = [int(v) for v in requested_auth if int(v) in allowed_auth_all]
        missing = [int(v) for v in requested_auth if int(v) not in allowed_auth_all]
        if missing:
            print(
                f"Warning: skipping requested WT trace auth_resid(s) not present in WT contact set: {missing}"
            )
        if not auth_resids:
            raise ValueError(
                "None of the requested WT trace auth_resid(s) are present in the WT contact set."
            )
    else:
        auth_resids = [int(v) for v in wt_occ_trace["auth_resid"].tolist()]

    label_by_auth = {int(a): str(l) for a, l in zip(wt_occ_all["auth_resid"], wt_occ_all["label"])}
    mean_traces = _extract_mutation_mean_traces_mdtraj(
        metas=metas,
        mutation=wt,
        auth_resids=auth_resids,
        ligand_resname=str(ligand_resname),
        resid_offset=int(resid_offset),
        window_ns=float(window_ns),
        n_grid=int(n_grid),
    )

    xpos = np.arange(len(wt_occ_all), dtype=float)
    bar_colors = wt_occ_all["region_color"].tolist()
    color_by_auth = {int(a): c for a, c in zip(wt_occ_all["auth_resid"], bar_colors)}

    traces: list[tuple[int, str, np.ndarray, np.ndarray, np.ndarray]] = []
    for auth_resid in auth_resids:
        gx, gm, gs = mean_traces[int(auth_resid)]
        traces.append((int(auth_resid), str(label_by_auth[int(auth_resid)]), gx, gm, gs))

    fig, (ax_top, ax_bot) = plt.subplots(
        2,
        1,
        figsize=(14.0, 9.0),
        gridspec_kw={"height_ratios": [1.35, 1.0]},
        constrained_layout=True,
    )

    for auth_resid, label, gx, gm, gs in traces:
        color = color_by_auth.get(int(auth_resid), "#1f77b4")
        ax_top.plot(
            gx,
            gm,
            color=color,
            linewidth=2.1,
            alpha=0.95,
            label=label,
        )
        lo = gm - gs
        hi = gm + gs
        ok = np.isfinite(lo) & np.isfinite(hi)
        if np.any(ok):
            ax_top.fill_between(gx[ok], lo[ok], hi[ok], color=color, alpha=0.16, linewidth=0)

    ax_top.axhline(
        float(contact_cutoff),
        color="#666666",
        linestyle=":",
        linewidth=1.2,
        label=f"{contact_cutoff:.1f} Å contact cutoff",
    )
    ax_top.set_xlim(0.0, float(window_ns))
    ax_top.set_xlabel("Time (ns)", fontsize=18)
    ax_top.set_ylabel("Min residue-DOR\n distance (Å)", fontsize=22)
    ax_top.tick_params(axis="both", labelsize=18)
    ax_top.spines["top"].set_visible(False)
    ax_top.spines["right"].set_visible(False)
    ax_top.grid(alpha=0.22, linestyle=":")
    ax_top.legend(loc="upper right", frameon=True, fontsize=15)

    means = wt_occ_all["occupancy_mean"].to_numpy(dtype=float)
    sems = wt_occ_all["occupancy_sem"].to_numpy(dtype=float)
    ax_bot.bar(
        xpos,
        means,
        yerr=sems,
        color=bar_colors,
        edgecolor="#333333",
        linewidth=0.4,
        capsize=3.0,
        alpha=0.95,
    )
    ax_bot.set_xticks(xpos, wt_occ_all["label"].tolist(), rotation=45, ha="right")
    ax_bot.tick_params(axis="both", labelsize=18)
    ax_bot.spines["top"].set_visible(False)
    ax_bot.spines["right"].set_visible(False)
    ax_bot.set_ylim(0.0, 1.05)
    ax_bot.set_ylabel("Mean occupancy", fontsize=22)
    ax_bot.set_xlabel("")
    ax_bot.grid(axis="y", alpha=0.2, linestyle=":")

    region_spans = (
        wt_occ_all.assign(xpos=xpos)
        .groupby(["region_order", "pocket_region"], as_index=False)
        .agg(xmin=("xpos", "min"), xmax=("xpos", "max"))
        .sort_values(["region_order", "xmin"])
    )
    for _, row in region_spans.iterrows():
        xmin = float(row["xmin"])
        xmax = float(row["xmax"])
        xmid = 0.5 * (xmin + xmax)
        span_width = max(1, int(round(xmax - xmin + 1.0)))
        region_label = "\n".join(
            textwrap.wrap(
                str(row["pocket_region"]),
                width=max(4, span_width * 4),
                break_long_words=False,
            )
        )
        ax_bot.text(
            xmid,
            1.02,
            region_label,
            transform=ax_bot.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=11,
            clip_on=False,
        )
        ax_bot.axvspan(xmin - 0.5, xmax + 0.5, color="#000000", alpha=0.025, zorder=0)
    for edge in region_spans["xmax"].tolist()[:-1]:
        ax_bot.axvline(float(edge) + 0.5, color="#666666", linestyle="--", linewidth=0.8, alpha=0.6)

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=320, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output_png}")


def _run_single_triplet(
    trip: tuple[str, str, str],
    metas: list[ReplicateMeta],
    ligand_resname: str,
    resid_offset: int,
    contact_cutoff: float,
    window_ns: float,
    trace_grid_points: int,
    min_wt_pooled_occ: float,
    min_neg_pooled_occ: float,
    max_can_pooled_occ: float,
    min_any_mean_occ_display: float,
    fold_map: dict[str, float],
    plots_dir: Path,
    output_prefix: str,
) -> dict[str, object]:
    rep_occ, mut_occ, timing_df = _compute_triplet_contact_stats(
        metas=metas,
        mutation_triplet=trip,
        ligand_resname=str(ligand_resname),
        resid_offset=int(resid_offset),
        contact_cutoff=float(contact_cutoff),
        window_ns=float(window_ns),
    )
    rep_occ["triplet"] = "|".join(trip)
    mut_occ["triplet"] = "|".join(trip)
    timing_df["triplet"] = "|".join(trip)

    story = _choose_story_residue(
        mut_occ=mut_occ,
        triplet=trip,
        min_wt_pooled=float(min_wt_pooled_occ),
        min_neg_pooled=float(min_neg_pooled_occ),
        max_can_pooled=float(max_can_pooled_occ),
        min_any_mean_occ=float(min_any_mean_occ_display),
    )
    auth_resid = int(story["auth_resid"])
    traj_resid = int(story["traj_resid"])
    story_res_map, story_has_pos_mut = _infer_resname_map_for_position(
        triplet=trip,
        auth_resid=int(auth_resid),
    )

    mean_traces: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    trace_rows: list[dict[str, object]] = []
    trace_rep_audit_rows: list[dict[str, object]] = []
    for mut in trip:
        gx, gm, gs, rep_audit = _extract_mutation_mean_trace(
            metas=metas,
            mutation=mut,
            auth_resid=int(auth_resid),
            ligand_resname=str(ligand_resname),
            resid_offset=int(resid_offset),
            window_ns=float(window_ns),
            n_grid=int(trace_grid_points),
        )
        mean_traces[mut] = (gx, gm, gs)
        for r in rep_audit:
            trace_rep_audit_rows.append({"triplet": "|".join(trip), "auth_resid": int(auth_resid), **r})
        for xi, yi, si in zip(gx, gm, gs):
            trace_rows.append(
                {
                    "triplet": "|".join(trip),
                    "mutation": mut,
                    "auth_resid": int(auth_resid),
                    "traj_resid": int(traj_resid),
                    "time_ns": float(xi),
                    "distance_mean_angstrom": float(yi),
                    "distance_sem_angstrom": float(si),
                    "dor_fold": float(fold_map.get(_normalize_mutation_token(mut), np.nan)),
                }
            )

    tag = f"{trip[0]}_{trip[1]}_{trip[2]}".replace("+", "_")
    out_png = plots_dir / f"{output_prefix}_{tag}.png"
    with _PLOT_LOCK:
        _plot_triplet_figure(
            triplet=trip,
            mut_occ=mut_occ,
            story_row=story,
            mean_traces=mean_traces,
            fold_map=fold_map,
            contact_cutoff=float(contact_cutoff),
            min_any_mean_occ_display=float(min_any_mean_occ_display),
            output_png=out_png,
        )

    summary_row = {
        "triplet": "|".join(trip),
        "wt": trip[0],
        "negative": trip[1],
        "canonical": trip[2],
        "wt_dor_fold": float(fold_map.get(_normalize_mutation_token(trip[0]), np.nan)),
        "negative_dor_fold": float(fold_map.get(_normalize_mutation_token(trip[1]), np.nan)),
        "canonical_dor_fold": float(fold_map.get(_normalize_mutation_token(trip[2]), np.nan)),
        "story_traj_resid": int(traj_resid),
        "story_auth_resid": int(auth_resid),
        "story_resname": str(story_res_map.get(trip[0], "UNK")),
        "story_label": _position_label(int(auth_resid)),
        "story_position_has_mutation": bool(story_has_pos_mut),
        "story_resname_wt": str(story_res_map.get(trip[0], "UNK")),
        "story_resname_negative": str(story_res_map.get(trip[1], "UNK")),
        "story_resname_canonical": str(story_res_map.get(trip[2], "UNK")),
        "story_direction": str(story.get("direction", "")),
        "story_score": float(pd.to_numeric(pd.Series([story.get("score")]), errors="coerce").iloc[0]),
        "story_occ_wt_pooled": float(pd.to_numeric(pd.Series([story.get(trip[0])]), errors="coerce").iloc[0]),
        "story_occ_negative_pooled": float(pd.to_numeric(pd.Series([story.get(trip[1])]), errors="coerce").iloc[0]),
        "story_occ_canonical_pooled": float(pd.to_numeric(pd.Series([story.get(trip[2])]), errors="coerce").iloc[0]),
        "output_png": str(out_png),
    }
    return {
        "triplet": trip,
        "rep_occ": rep_occ,
        "mut_occ": mut_occ,
        "timing_df": timing_df,
        "trace_rows": trace_rows,
        "trace_rep_audit_rows": trace_rep_audit_rows,
        "summary_row": summary_row,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate pooled-contact + mean-trace story figures for mutation triplets.")
    parser.add_argument("--manifest", type=Path, default=Path("manifests/md_manifest.csv"))
    parser.add_argument("--susceptibility-xlsx", type=Path, default=Path("data/DRM-susceptibilities.csv.xlsx"))
    parser.add_argument("--ligand-resname", type=str, default="2KW")
    parser.add_argument("--resid-offset", type=int, default=-3, help="auth = traj - resid_offset")
    parser.add_argument("--contact-cutoff", type=float, default=4.0)
    parser.add_argument("--window-ns", type=float, default=100.0)
    parser.add_argument("--trace-grid-points", type=int, default=500)
    parser.add_argument("--min-wt-pooled-occ", type=float, default=0.45)
    parser.add_argument("--min-neg-pooled-occ", type=float, default=0.45)
    parser.add_argument("--max-can-pooled-occ", type=float, default=0.70)
    parser.add_argument(
        "--min-any-mean-occ-display",
        type=float,
        default=0.5,
        help="Display/select residues only when contact_mean is strictly greater than this value in at least one genotype.",
    )
    parser.add_argument(
        "--triplets",
        type=str,
        default=(
            "WT,V106M,V106A;"
            "WT,G190A,G190E;"
            "WT,K103N,K103N+P225H;"
            "WT,K103N,K103N+M230L;"
            "WT,V106A,V106A+P225H;"
            "WT,V106I,V106I+F227C;"
            "WT,V106I,V106A+P225H;"
            "WT,V106I,V106A+L234I;"
            "WT,V106I,V106A;"
            "WT,V106A,V106A+L234I;"
            "WT,F227C,A98G+F227C;"
            "WT,F227C,V106I+F227C;"
            "WT,Y181C,Y188L;"
            "WT,K103N,Y318F"
        ),
        help="Semicolon-separated triplets, each as WT,NEG,CAN.",
    )
    parser.add_argument("--output-base-dir", type=Path, default=Path("results/analysis/triplet_contact_story_100ns"))
    parser.add_argument("--output-prefix", type=str, default="triplet_story_100ns")
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Number of parallel worker threads for per-triplet generation (1 = serial).",
    )
    parser.add_argument(
        "--wt-only-mutation",
        type=str,
        default="",
        help="Optional mutation name for a WT-only contacts figure computed from that mutation only.",
    )
    parser.add_argument(
        "--wt-only-only",
        action="store_true",
        help="If set, generate only the WT-only contacts figure and skip triplet story outputs.",
    )
    parser.add_argument(
        "--wt-only-triplet",
        type=str,
        default="",
        help="Optional single triplet (WT,NEG,CAN) for an additional WT-only contacts figure.",
    )
    parser.add_argument(
        "--wt-only-output-name",
        type=str,
        default="100_ns_WT_contacts.png",
        help="Output filename (inside output-base-dir/plots) for WT-only contacts figure.",
    )
    parser.add_argument(
        "--wt-only-min-mean-occ",
        type=float,
        default=0.5,
        help="Minimum WT contact_mean required for inclusion in WT-only figure.",
    )
    parser.add_argument(
        "--wt-only-trace-auth-resids",
        type=str,
        default="105,179",
        help="Comma-separated auth residue numbers to show in WT-only top-panel traces.",
    )
    args = parser.parse_args()

    triplets: list[tuple[str, str, str]] = []
    if not bool(args.wt_only_only):
        triplets = _parse_triplets(args.triplets)
    needed_mutations = {m for t in triplets for m in t}
    wt_only_triplet: tuple[str, str, str] | None = None
    wt_only_mutation = str(args.wt_only_mutation).strip()
    if wt_only_mutation:
        needed_mutations.add(wt_only_mutation)
    if str(args.wt_only_triplet).strip():
        wt_only_parsed = _parse_triplets(str(args.wt_only_triplet))
        if len(wt_only_parsed) != 1:
            raise ValueError("--wt-only-triplet must contain exactly one triplet.")
        wt_only_triplet = wt_only_parsed[0]
        needed_mutations.update(set(wt_only_triplet))
    if bool(args.wt_only_only) and not wt_only_mutation:
        raise ValueError("--wt-only-only requires --wt-only-mutation.")
    wt_only_trace_auth_resids = _parse_int_csv(str(args.wt_only_trace_auth_resids))

    fold_map = _load_dor_fold_map(args.susceptibility_xlsx)
    metas = _load_replicate_meta(args.manifest, needed_mutations=needed_mutations)
    if not metas:
        raise ValueError("No usable replicate metadata found.")

    output_base = args.output_base_dir
    plots_dir = output_base / "plots"
    tables_dir = output_base / "tables"
    config_dir = output_base / "config"
    for d in [plots_dir, tables_dir, config_dir]:
        d.mkdir(parents=True, exist_ok=True)

    if wt_only_mutation:
        wt_triplet = (wt_only_mutation, wt_only_mutation, wt_only_mutation)
        _, wt_mut_occ, _ = _compute_triplet_contact_stats(
            metas=metas,
            mutation_triplet=wt_triplet,
            ligand_resname=str(args.ligand_resname),
            resid_offset=int(args.resid_offset),
            contact_cutoff=float(args.contact_cutoff),
            window_ns=float(args.window_ns),
        )
        wt_only_png = plots_dir / str(args.wt_only_output_name)
        _plot_wt_contacts_figure(
            wt_mutation=wt_only_mutation,
            mut_occ=wt_mut_occ,
            metas=metas,
            ligand_resname=str(args.ligand_resname),
            resid_offset=int(args.resid_offset),
            window_ns=float(args.window_ns),
            n_grid=int(args.trace_grid_points),
            contact_cutoff=float(args.contact_cutoff),
            min_wt_mean_occ=float(args.wt_only_min_mean_occ),
            wt_trace_auth_resids=wt_only_trace_auth_resids,
            output_png=wt_only_png,
        )
        if bool(args.wt_only_only):
            return 0

    summary_rows: list[dict[str, object]] = []
    all_mut_occ: list[pd.DataFrame] = []
    all_rep_occ: list[pd.DataFrame] = []
    all_timing: list[pd.DataFrame] = []
    all_trace_rows: list[dict[str, object]] = []
    trace_rep_audit_rows: list[dict[str, object]] = []

    triplet_results: list[dict[str, object]] = []
    workers = max(1, int(args.workers))
    if triplets:
        if workers <= 1 or len(triplets) <= 1:
            for trip in triplets:
                triplet_results.append(
                    _run_single_triplet(
                        trip=trip,
                        metas=metas,
                        ligand_resname=str(args.ligand_resname),
                        resid_offset=int(args.resid_offset),
                        contact_cutoff=float(args.contact_cutoff),
                        window_ns=float(args.window_ns),
                        trace_grid_points=int(args.trace_grid_points),
                        min_wt_pooled_occ=float(args.min_wt_pooled_occ),
                        min_neg_pooled_occ=float(args.min_neg_pooled_occ),
                        max_can_pooled_occ=float(args.max_can_pooled_occ),
                        min_any_mean_occ_display=float(args.min_any_mean_occ_display),
                        fold_map=fold_map,
                        plots_dir=plots_dir,
                        output_prefix=str(args.output_prefix),
                    )
                )
        else:
            max_workers = min(len(triplets), workers)
            future_map = {}
            with ThreadPoolExecutor(max_workers=max_workers) as ex:
                for trip in triplets:
                    fut = ex.submit(
                        _run_single_triplet,
                        trip,
                        metas,
                        str(args.ligand_resname),
                        int(args.resid_offset),
                        float(args.contact_cutoff),
                        float(args.window_ns),
                        int(args.trace_grid_points),
                        float(args.min_wt_pooled_occ),
                        float(args.min_neg_pooled_occ),
                        float(args.max_can_pooled_occ),
                        float(args.min_any_mean_occ_display),
                        fold_map,
                        plots_dir,
                        str(args.output_prefix),
                    )
                    future_map[fut] = trip
                triplet_result_by_trip: dict[tuple[str, str, str], dict[str, object]] = {}
                for fut in as_completed(future_map):
                    trip = future_map[fut]
                    triplet_result_by_trip[trip] = fut.result()
                triplet_results = [triplet_result_by_trip[t] for t in triplets]

    for res in triplet_results:
        all_rep_occ.append(res["rep_occ"])
        all_mut_occ.append(res["mut_occ"])
        all_timing.append(res["timing_df"])
        all_trace_rows.extend(res["trace_rows"])
        trace_rep_audit_rows.extend(res["trace_rep_audit_rows"])
        summary_rows.append(res["summary_row"])

    if wt_only_triplet is not None and not bool(args.wt_only_only):
        wt_target = tuple(wt_only_triplet)
        match = None
        for res in triplet_results:
            if tuple(res["triplet"]) == wt_target:
                match = res
                break
        if match is None:
            raise ValueError(f"WT-only triplet not found in generated triplets: {wt_target}")
        wt_only_png = plots_dir / str(args.wt_only_output_name)
        _plot_wt_contacts_figure(
            wt_mutation=wt_target[0],
            mut_occ=match["mut_occ"],
            metas=metas,
            ligand_resname=str(args.ligand_resname),
            resid_offset=int(args.resid_offset),
            window_ns=float(args.window_ns),
            n_grid=int(args.trace_grid_points),
            contact_cutoff=float(args.contact_cutoff),
            min_wt_mean_occ=float(args.wt_only_min_mean_occ),
            wt_trace_auth_resids=wt_only_trace_auth_resids,
            output_png=wt_only_png,
        )

    summary_df = pd.DataFrame(summary_rows)
    rep_df = pd.concat(all_rep_occ, ignore_index=True) if all_rep_occ else pd.DataFrame()
    mut_df = pd.concat(all_mut_occ, ignore_index=True) if all_mut_occ else pd.DataFrame()
    timing_df = pd.concat(all_timing, ignore_index=True) if all_timing else pd.DataFrame()
    trace_df = pd.DataFrame(all_trace_rows)
    trace_rep_audit_df = pd.DataFrame(trace_rep_audit_rows)

    summary_df = summary_df.rename(
        columns={
            "story_occ_wt_pooled": "story_contact_wt_pooled",
            "story_occ_negative_pooled": "story_contact_negative_pooled",
            "story_occ_canonical_pooled": "story_contact_canonical_pooled",
        }
    )
    rep_df = rep_df.rename(columns=lambda c: str(c).replace("occupancy", "contact"))
    mut_df = mut_df.rename(columns=lambda c: str(c).replace("occupancy", "contact"))

    summary_csv = tables_dir / "selection_summary.csv"
    rep_csv = tables_dir / "replicate_contact.csv"
    mut_csv = tables_dir / "mutation_contact.csv"
    trace_csv = tables_dir / "mean_traces.csv"
    trace_rep_csv = tables_dir / "trace_replicate_audit.csv"
    timing_csv = tables_dir / "timing_audit.csv"
    fold_csv = tables_dir / "fold_lookup.csv"

    summary_df.to_csv(summary_csv, index=False)
    rep_df.to_csv(rep_csv, index=False)
    mut_df.to_csv(mut_csv, index=False)
    trace_df.to_csv(trace_csv, index=False)
    trace_rep_audit_df.to_csv(trace_rep_csv, index=False)
    timing_df.to_csv(timing_csv, index=False)
    pd.DataFrame(sorted(fold_map.items()), columns=["mutation", "dor_fold"]).to_csv(fold_csv, index=False)

    (config_dir / "triplets.txt").write_text("\n".join(["|".join(t) for t in triplets]) + "\n")

    print(f"Saved {summary_csv}")
    print(f"Saved {rep_csv}")
    print(f"Saved {mut_csv}")
    print(f"Saved {trace_csv}")
    print(f"Saved {trace_rep_csv}")
    print(f"Saved {timing_csv}")
    print(f"Saved {fold_csv}")
    print(f"Saved {config_dir / 'triplets.txt'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
