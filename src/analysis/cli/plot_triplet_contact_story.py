#!/usr/bin/env python3
"""Generate 100 ns triplet contact-story figures using pooled occupancy and mean traces.

For each (WT, comparator, DRM) triplet this script:
1) Computes residue-level occupancy across all frames and all replicates (pooled).
2) Picks a story residue where canonical occupancy is lower than WT/comparator.
3) Plots mean distance traces (average over all replicates) for the story residue.
4) Adds a raw occupancy heatmap (pooled occupancy) for all contacted residues.
5) Labels legend with DOR fold-change values from DRM-susceptibilities sheet.
"""
from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


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

        if np.isfinite(ns_state) and ns_state > 0:
            total_ns = float(ns_state)
            timing_source = "state_csv"
        elif np.isfinite(ns_json) and ns_json > 0:
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


def _select_residue_sidechain(universe, traj_resid: int):
    sc = universe.select_atoms(f"protein and resid {int(traj_resid)} and not name H* and not backbone")
    if sc.n_atoms > 0:
        return sc
    gly = universe.select_atoms(f"protein and resid {int(traj_resid)} and name CA")
    if gly.n_atoms > 0:
        return gly
    return universe.select_atoms(f"protein and resid {int(traj_resid)} and not name H*")


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
        prot = u.select_atoms("(protein and not name H* and not backbone) or (protein and resname GLY and name CA)")
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
) -> pd.Series:
    wt, neg, can = triplet
    piv = (
        mut_occ.pivot_table(
            index=["traj_resid", "auth_resid", "resname"],
            columns="mutation",
            values="occupancy_pooled",
            aggfunc="mean",
        )
        .reset_index()
    )
    for c in [wt, neg, can]:
        if c not in piv.columns:
            piv[c] = np.nan
    piv["score"] = np.minimum(pd.to_numeric(piv[wt], errors="coerce"), pd.to_numeric(piv[neg], errors="coerce")) - pd.to_numeric(piv[can], errors="coerce")
    piv = piv.replace([np.inf, -np.inf], np.nan).dropna(subset=[wt, neg, can, "score"]).copy()
    if piv.empty:
        raise ValueError(f"No comparable residues across triplet {triplet}")

    constrained = piv[
        (pd.to_numeric(piv[wt], errors="coerce") >= float(min_wt_pooled))
        & (pd.to_numeric(piv[neg], errors="coerce") >= float(min_neg_pooled))
        & (pd.to_numeric(piv[can], errors="coerce") <= float(max_can_pooled))
        & (pd.to_numeric(piv["score"], errors="coerce") > 0.0)
    ].copy()
    if not constrained.empty:
        return constrained.sort_values(["score", wt, neg], ascending=[False, False, False]).iloc[0]
    return piv.sort_values("score", ascending=False).iloc[0]


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
        x, y = _extract_distance_trace(
            meta=m,
            auth_resid=int(auth_resid),
            ligand_resname=ligand_resname,
            resid_offset=int(resid_offset),
            window_ns=float(window_ns),
        )
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


def _format_fold_label(mutation: str, fold_map: dict[str, float]) -> str:
    key = _normalize_mutation_token(mutation)
    v = fold_map.get(key, np.nan)
    if np.isfinite(v):
        if abs(v - round(v)) < 1e-8:
            return f"{mutation} ({int(round(v))}x)"
        return f"{mutation} ({v:.1f}x)"
    return mutation


def _plot_triplet_figure(
    triplet: tuple[str, str, str],
    mut_occ: pd.DataFrame,
    story_row: pd.Series,
    mean_traces: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]],
    fold_map: dict[str, float],
    contact_cutoff: float,
    output_png: Path,
) -> None:
    import matplotlib.pyplot as plt

    wt, neg, can = triplet
    mut_order = [wt, neg, can]
    colors = {wt: "#4d4d4d", neg: "#1f77b4", can: "#d62728"}

    all_res = (
        mut_occ[["traj_resid", "auth_resid", "resname"]]
        .drop_duplicates()
        .sort_values(["auth_resid", "traj_resid"])
        .copy()
    )
    all_res["label"] = all_res.apply(lambda r: _residue_label(int(r["auth_resid"]), str(r["resname"])), axis=1)
    col_order = all_res["label"].tolist()
    key_lookup = all_res.copy()
    key_lookup["key"] = list(zip(key_lookup["traj_resid"], key_lookup["auth_resid"], key_lookup["resname"]))
    map_lab = dict(zip(key_lookup["key"], key_lookup["label"]))
    tmp = mut_occ.copy()
    tmp["label"] = [map_lab[(int(a), int(b), str(c))] for a, b, c in zip(tmp["traj_resid"], tmp["auth_resid"], tmp["resname"])]
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

    ax_top.axhline(float(contact_cutoff), color="#666666", linestyle=":", linewidth=1.1, label=f"{contact_cutoff:.1f} A cutoff")
    ax_top.set_xlim(0.0, xmax if xmax > 0 else 100.0)
    ax_top.set_xlabel("Time (ns)")
    ax_top.set_ylabel("Min sidechain-DOR distance (A)")
    story_label = _residue_label(int(story_row["auth_resid"]), str(story_row.get("resname", "")))
    ax_top.set_title(f"Selected Story Residue: {story_label} (canonical-loss pattern)")
    ax_top.grid(alpha=0.22, linestyle=":")
    ax_top.legend(loc="upper right", frameon=True, fontsize=8)

    score = float(pd.to_numeric(pd.Series([story_row.get("score")]), errors="coerce").iloc[0])
    wt_occ = float(pd.to_numeric(pd.Series([story_row.get(wt)]), errors="coerce").iloc[0])
    neg_occ = float(pd.to_numeric(pd.Series([story_row.get(neg)]), errors="coerce").iloc[0])
    can_occ = float(pd.to_numeric(pd.Series([story_row.get(can)]), errors="coerce").iloc[0])
    txt = (
        f"Pooled occupancy at {story_label}: {wt}={wt_occ:.2f}, {neg}={neg_occ:.2f}, {can}={can_occ:.2f}\n"
        f"score=min(WT,comp)-can={score:.2f}"
    )
    ax_top.text(0.01, 0.98, txt, transform=ax_top.transAxes, va="top", ha="left", fontsize=8.5, bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "#cccccc"})

    arr = hm.to_numpy(dtype=float)
    im = ax_bot.imshow(arr, aspect="auto", cmap="viridis", vmin=0.0, vmax=1.0)
    ax_bot.set_yticks(np.arange(len(mut_order)), [_format_fold_label(m, fold_map=fold_map) for m in mut_order])
    ax_bot.set_xticks(np.arange(len(col_order)), col_order, rotation=45, ha="right")
    ax_bot.set_xlabel("Residue")
    ax_bot.set_ylabel("Mutation")
    ax_bot.set_title("Raw Occupancy Heatmap (first 100 ns, all replicates, all frames; pooled)")
    cbar = fig.colorbar(im, ax=ax_bot, fraction=0.046, pad=0.02)
    cbar.set_label("Occupancy")

    fig.suptitle(f"Triplet Contact Story: {_format_fold_label(wt, fold_map)} vs {_format_fold_label(neg, fold_map)} vs {_format_fold_label(can, fold_map)}", fontsize=13.5, fontweight="bold")
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=320, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output_png}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate pooled-occupancy + mean-trace story figures for mutation triplets.")
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
        "--triplets",
        type=str,
        default=(
            "WT,V106M,V106A;"
            "WT,G190A,G190E;"
            "WT,K103N,K103N+P225H;"
            "WT,K103N,K103N+M230L;"
            "WT,V106A,V106A+P225H;"
            "WT,V106I,V106I+F227C;"
            "WT,V106A,V106A+L234I;"
            "WT,F227C,A98G+F227C"
        ),
        help="Semicolon-separated triplets, each as WT,NEG,CAN.",
    )
    parser.add_argument("--output-base-dir", type=Path, default=Path("results/analysis/triplet_contact_story_100ns"))
    parser.add_argument("--output-prefix", type=str, default="triplet_story_100ns")
    args = parser.parse_args()

    triplets = _parse_triplets(args.triplets)
    needed_mutations = {m for t in triplets for m in t}

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

    summary_rows: list[dict[str, object]] = []
    all_mut_occ: list[pd.DataFrame] = []
    all_rep_occ: list[pd.DataFrame] = []
    all_timing: list[pd.DataFrame] = []
    all_trace_rows: list[dict[str, object]] = []
    trace_rep_audit_rows: list[dict[str, object]] = []

    for trip in triplets:
        rep_occ, mut_occ, timing_df = _compute_triplet_contact_stats(
            metas=metas,
            mutation_triplet=trip,
            ligand_resname=str(args.ligand_resname),
            resid_offset=int(args.resid_offset),
            contact_cutoff=float(args.contact_cutoff),
            window_ns=float(args.window_ns),
        )
        rep_occ["triplet"] = "|".join(trip)
        mut_occ["triplet"] = "|".join(trip)
        timing_df["triplet"] = "|".join(trip)
        all_rep_occ.append(rep_occ)
        all_mut_occ.append(mut_occ)
        all_timing.append(timing_df)

        story = _choose_story_residue(
            mut_occ=mut_occ,
            triplet=trip,
            min_wt_pooled=float(args.min_wt_pooled_occ),
            min_neg_pooled=float(args.min_neg_pooled_occ),
            max_can_pooled=float(args.max_can_pooled_occ),
        )
        auth_resid = int(story["auth_resid"])
        traj_resid = int(story["traj_resid"])

        mean_traces: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
        for mut in trip:
            gx, gm, gs, rep_audit = _extract_mutation_mean_trace(
                metas=metas,
                mutation=mut,
                auth_resid=int(auth_resid),
                ligand_resname=str(args.ligand_resname),
                resid_offset=int(args.resid_offset),
                window_ns=float(args.window_ns),
                n_grid=int(args.trace_grid_points),
            )
            mean_traces[mut] = (gx, gm, gs)
            for r in rep_audit:
                trace_rep_audit_rows.append({"triplet": "|".join(trip), "auth_resid": int(auth_resid), **r})
            for xi, yi, si in zip(gx, gm, gs):
                all_trace_rows.append(
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
        out_png = plots_dir / f"{args.output_prefix}_{tag}.png"
        _plot_triplet_figure(
            triplet=trip,
            mut_occ=mut_occ,
            story_row=story,
            mean_traces=mean_traces,
            fold_map=fold_map,
            contact_cutoff=float(args.contact_cutoff),
            output_png=out_png,
        )

        summary_rows.append(
            {
                "triplet": "|".join(trip),
                "wt": trip[0],
                "negative": trip[1],
                "canonical": trip[2],
                "wt_dor_fold": float(fold_map.get(_normalize_mutation_token(trip[0]), np.nan)),
                "negative_dor_fold": float(fold_map.get(_normalize_mutation_token(trip[1]), np.nan)),
                "canonical_dor_fold": float(fold_map.get(_normalize_mutation_token(trip[2]), np.nan)),
                "story_traj_resid": int(traj_resid),
                "story_auth_resid": int(auth_resid),
                "story_resname": str(story.get("resname", "")),
                "story_label": _residue_label(int(auth_resid), str(story.get("resname", ""))),
                "story_score": float(pd.to_numeric(pd.Series([story.get("score")]), errors="coerce").iloc[0]),
                "story_occ_wt_pooled": float(pd.to_numeric(pd.Series([story.get(trip[0])]), errors="coerce").iloc[0]),
                "story_occ_negative_pooled": float(pd.to_numeric(pd.Series([story.get(trip[1])]), errors="coerce").iloc[0]),
                "story_occ_canonical_pooled": float(pd.to_numeric(pd.Series([story.get(trip[2])]), errors="coerce").iloc[0]),
                "output_png": str(out_png),
            }
        )

    summary_df = pd.DataFrame(summary_rows)
    rep_df = pd.concat(all_rep_occ, ignore_index=True) if all_rep_occ else pd.DataFrame()
    mut_df = pd.concat(all_mut_occ, ignore_index=True) if all_mut_occ else pd.DataFrame()
    timing_df = pd.concat(all_timing, ignore_index=True) if all_timing else pd.DataFrame()
    trace_df = pd.DataFrame(all_trace_rows)
    trace_rep_audit_df = pd.DataFrame(trace_rep_audit_rows)

    summary_csv = tables_dir / "selection_summary.csv"
    rep_csv = tables_dir / "replicate_occupancy.csv"
    mut_csv = tables_dir / "mutation_occupancy.csv"
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
