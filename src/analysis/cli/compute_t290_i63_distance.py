#!/usr/bin/env python3
"""Compute T290-I63 C-alpha distance across apo+holo trajectories.

Outputs:
  - results/reference_t290_i63_distance.csv
  - results/t290_i63_distance_dynamics.csv
  - results/t290_i63_distance_summary.csv
  - results/apo_t290_i63_distance.csv
  - results/plots/apo_t290_i63_distance_{timeseries,distribution,vs_fold}.png

"""

from __future__ import annotations

import json
import re
from collections import Counter
from pathlib import Path

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser

REPO = Path(__file__).resolve().parents[3]
RESULTS = REPO / "results"
PLOTS = RESULTS / "plots"

HOLO_MANIFEST = RESULTS / "md_manifest.csv"
APO_MANIFEST = RESULTS / "apo_md_manifest.csv"
SUSCEPTIBILITY_XLSX = REPO / "data" / "DRM-susceptibilities.csv.xlsx"

OUT_REF = RESULTS / "reference_t290_i63_distance.csv"
OUT_DYNAMICS = RESULTS / "t290_i63_distance_dynamics.csv"
OUT_SUMMARY = RESULTS / "t290_i63_distance_summary.csv"
OUT_APO = RESULTS / "apo_t290_i63_distance.csv"

RESID_OFFSET = -3
I63_CANON = 63
T290_CANON = 290
ATOM_NAME = "CA"


def _normalize_mutation_label(text: str) -> str:
    t = str(text or "").strip()
    if not t:
        return ""
    t = re.sub(r"\s+", "", t)
    t = t.replace(",", "+")
    t = re.sub(r"\++", "+", t)
    return t


def _load_fold_map() -> dict[str, float]:
    fold: dict[str, float] = {}
    # Baseline from manifests (fallback values).
    for mf in (HOLO_MANIFEST, APO_MANIFEST):
        if not mf.exists():
            continue
        df = pd.read_csv(mf)
        for _, row in df.drop_duplicates("mutation").iterrows():
            mut = _normalize_mutation_label(str(row.get("mutation", "")))
            if not mut:
                continue
            val = row.get("fold_reduction")
            if pd.notna(val):
                try:
                    fold[mut] = float(val)
                except Exception:
                    pass
    # Override with source-of-truth DOR values from susceptibility table.
    if SUSCEPTIBILITY_XLSX.exists():
        sdf = pd.read_excel(SUSCEPTIBILITY_XLSX)
        if len(sdf.columns) >= 3:
            sdf = sdf.iloc[:, :3].copy()
            sdf.columns = ["mutation", "rpv_fold", "dor_fold"]
            for _, row in sdf.iterrows():
                mut = _normalize_mutation_label(row.get("mutation"))
                if not mut or mut.lower() == "mutations":
                    continue
                v = row.get("dor_fold")
                if pd.notna(v):
                    try:
                        fold[mut] = float(v)
                    except Exception:
                        pass
    # Keep WT explicit for legends.
    fold.setdefault("WT", 1.0)
    return fold


def _format_fold(fold: float | int | None) -> str:
    try:
        f = float(fold)  # type: ignore[arg-type]
    except Exception:
        return "?"
    if not np.isfinite(f):
        return "?"
    if abs(f - round(f)) < 1e-6:
        return f"{int(round(f))}x"
    return f"{f:.1f}x"


def _remap_to_local_workspace(candidate: Path | None) -> Path | None:
    if candidate is None:
        return None
    if candidate.exists():
        return candidate
    marker = "nnrti-mechanisms/"
    text = str(candidate)
    if marker not in text:
        return candidate
    rel = text.split(marker, 1)[1]
    mapped = REPO / rel
    if mapped.exists():
        return mapped
    return candidate


def _traj_paths_from_row(row: pd.Series) -> tuple[Path | None, Path | None]:
    out_json = Path(str(row.get("output_json", "")).strip())
    if not out_json.exists():
        return None, None
    data = json.loads(out_json.read_text())
    topo = Path(str(data.get("analysis_topology_pdb") or "").strip())
    dcd = Path(str(data.get("analysis_dcd") or "").strip())
    topo = _remap_to_local_workspace(topo)
    dcd = _remap_to_local_workspace(dcd)
    if topo is None or dcd is None:
        return None, None
    if not topo.exists() or not dcd.exists():
        return None, None
    return topo, dcd


def _infer_total_ns(output_json_path: Path) -> float:
    try:
        j = json.loads(output_json_path.read_text())
        steps = int(j.get("md_production_steps_completed") or j.get("md_production_steps") or 0)
        if steps > 0:
            return steps * 2.0 / 1_000_000.0
    except Exception:
        pass
    m = re.match(r"^(.+)_rep(\d{2}).*\.json$", output_json_path.name)
    if m:
        state_csv = output_json_path.parent / f"{m.group(1)}_rep{m.group(2)}_md_state.csv"
        if state_csv.exists():
            try:
                sdf = pd.read_csv(state_csv)
                for col in ('#"Step"', "Step"):
                    if col in sdf.columns:
                        steps = pd.to_numeric(sdf[col], errors="coerce").dropna()
                        if not steps.empty:
                            return float(steps.max()) * 2.0 / 1_000_000.0
            except Exception:
                pass
    return 100.0


def _p66_chain_for_cif(cif_path: Path) -> str:
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(cif_path.stem, str(cif_path))
    model = next(structure.get_models())
    chain_counts: dict[str, int] = {}
    for ch in model:
        n_ca = sum(1 for r in ch if r.id[0] == " " and "CA" in r)
        chain_counts[ch.id] = n_ca
    if not chain_counts:
        raise ValueError(f"No protein CA chains in {cif_path}")
    return max(chain_counts, key=chain_counts.get)


def _reference_distance(cif_path: Path, resid_i63: int, resid_t290: int, atom_name: str) -> float:
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(cif_path.stem, str(cif_path))
    model = next(structure.get_models())
    p66 = _p66_chain_for_cif(cif_path)
    chain = model[p66]
    ai = chain[(" ", resid_i63, " ")][atom_name].coord
    at = chain[(" ", resid_t290, " ")][atom_name].coord
    return float(np.linalg.norm(at - ai))


def _p66_segid(universe) -> str:
    prot_ca = universe.select_atoms("protein and name CA")
    if prot_ca.n_atoms == 0:
        return ""
    cnt = Counter(prot_ca.segids.tolist())
    return max(cnt, key=cnt.get)


def _process_replicate(row: pd.Series, leg: str, frame_stride: int = 1) -> list[dict]:
    import MDAnalysis as mda

    topo, dcd = _traj_paths_from_row(row)
    if topo is None or dcd is None:
        return []

    mutation = str(row["mutation"])
    safe_label = str(row.get("safe_label", mutation))
    replicate = int(row["replicate"])
    output_json = Path(str(row["output_json"]))
    total_ns = _infer_total_ns(output_json) if output_json.exists() else 100.0

    u = mda.Universe(str(topo), str(dcd))
    seg = _p66_segid(u)
    seg_filter = f" and segid {seg}" if seg else ""
    resid_i63 = I63_CANON + RESID_OFFSET
    resid_t290 = T290_CANON + RESID_OFFSET
    sel_i63 = u.select_atoms(
        f"protein{seg_filter} and resid {resid_i63} and name {ATOM_NAME}"
    )
    sel_t290 = u.select_atoms(
        f"protein{seg_filter} and resid {resid_t290} and name {ATOM_NAME}"
    )
    if sel_i63.n_atoms == 0 or sel_t290.n_atoms == 0:
        return []

    n_frames = len(u.trajectory)
    out: list[dict] = []
    for ts in u.trajectory[:: max(1, int(frame_stride))]:
        d = float(np.linalg.norm(sel_t290.positions[0] - sel_i63.positions[0]))
        time_ns = (float(ts.frame) / max(1, n_frames - 1)) * total_ns
        out.append(
            {
                "leg": leg,
                "mutation": mutation,
                "safe_label": safe_label,
                "replicate": replicate,
                "frame": int(ts.frame),
                "time_ns": time_ns,
                "t290_i63_ca_distance_angstrom": d,
            }
        )
    return out


def _plot_apo(
    apo_df: pd.DataFrame,
    fold_map: dict[str, float],
    ref_1dlo: float,
    ref_4ncg: float,
) -> None:
    if apo_df.empty:
        return

    PLOTS.mkdir(parents=True, exist_ok=True)
    mutations = sorted(
        apo_df["mutation"].unique(),
        key=lambda m: (m != "WT", fold_map.get(m, 9999.0), m),
    )
    colors = dict(zip(mutations, cm.RdYlGn_r(np.linspace(0.1, 0.9, len(mutations)))))

    # Midpoint threshold to classify 1DLO-like vs 4NCG-like.
    threshold = 0.5 * (ref_1dlo + ref_4ncg)
    one_dlo_is_lower = ref_1dlo < ref_4ncg

    # Timeseries
    fig, axes = plt.subplots(len(mutations), 1, figsize=(12, max(2.2 * len(mutations), 4)), sharex=False)
    if len(mutations) == 1:
        axes = [axes]
    for ax, mut in zip(axes, mutations):
        sub = apo_df[apo_df["mutation"] == mut]
        for rep, grp in sub.groupby("replicate"):
            g = grp.sort_values("time_ns")
            ax.plot(
                g["time_ns"],
                g["t290_i63_ca_distance_angstrom"],
                linewidth=0.8,
                alpha=0.75,
                color=colors[mut],
                label=f"rep{rep}",
            )
        ax.axhline(ref_1dlo, color="navy", linestyle="--", linewidth=1.0, alpha=0.9, label="1DLO (Apo)")
        ax.axhline(ref_4ncg, color="darkgreen", linestyle="--", linewidth=1.0, alpha=0.9, label="4NCG (Holo)")
        ax.axhline(threshold, color="gray", linestyle=":", linewidth=0.9, alpha=0.7, label=f"midpoint ({threshold:.2f} Å)")
        ax.set_ylabel("Distance (Å)", fontsize=8)
        fold = fold_map.get(mut, np.nan)
        fold_txt = _format_fold(fold)
        ax.set_title(f"{mut}  (fold={fold_txt})", fontsize=9, loc="left", fontweight="bold")
        ax.grid(alpha=0.2, linestyle=":")
        ax.legend(fontsize=6, frameon=False, ncol=4, loc="upper right")
    axes[-1].set_xlabel("Time (ns)")
    fig.suptitle("T290-I63 Cα distance — apo simulations", fontsize=11, fontweight="bold")
    fig.tight_layout()
    out_timeseries = PLOTS / "apo_t290_i63_distance_timeseries.png"
    fig.savefig(out_timeseries, dpi=180, bbox_inches="tight")
    plt.close(fig)

    # Distribution
    fig, ax = plt.subplots(figsize=(10, 4))
    vals_all = apo_df["t290_i63_ca_distance_angstrom"].to_numpy()
    bins = np.linspace(max(0.0, vals_all.min() - 1.0), vals_all.max() + 1.0, 80)
    for mut in mutations:
        vals = apo_df[apo_df["mutation"] == mut]["t290_i63_ca_distance_angstrom"].to_numpy()
        counts, edges = np.histogram(vals, bins=bins, density=True)
        centers = 0.5 * (edges[:-1] + edges[1:])
        fold = fold_map.get(mut, np.nan)
        fold_txt = _format_fold(fold)
        ax.plot(centers, counts, color=colors[mut], linewidth=1.5, label=f"{mut} ({fold_txt})")
        ax.fill_between(centers, counts, alpha=0.08, color=colors[mut])
    ax.axvline(ref_1dlo, color="navy", linestyle="--", linewidth=1.2, alpha=0.9, label="1DLO (Apo)")
    ax.axvline(ref_4ncg, color="darkgreen", linestyle="--", linewidth=1.2, alpha=0.9, label="4NCG (Holo)")
    ax.set_xlabel("T290-I63 Cα distance (Å)")
    ax.set_ylabel("Density")
    ax.set_title("T290-I63 distance distribution — apo simulations", fontweight="bold")
    ax.legend(fontsize=8, frameon=False)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    out_dist = PLOTS / "apo_t290_i63_distance_distribution.png"
    fig.savefig(out_dist, dpi=180, bbox_inches="tight")
    plt.close(fig)

    # Vs fold
    records: list[dict] = []
    for mut in mutations:
        sub = apo_df[apo_df["mutation"] == mut]
        for rep, grp in sub.groupby("replicate"):
            vals = grp["t290_i63_ca_distance_angstrom"].to_numpy()
            if one_dlo_is_lower:
                frac_1dlo_like = float(np.mean(vals < threshold))
            else:
                frac_1dlo_like = float(np.mean(vals > threshold))
            records.append(
                {
                    "mutation": mut,
                    "replicate": int(rep),
                    "fold": fold_map.get(mut, np.nan),
                    "mean_distance": float(np.mean(vals)),
                    "frac_1dlo_like": frac_1dlo_like,
                }
            )
    rdf = pd.DataFrame(records)

    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    for mut in mutations:
        sub = rdf[rdf["mutation"] == mut]
        if sub.empty:
            continue
        fold = sub["fold"].iloc[0]
        ax.scatter(sub["fold"], sub["frac_1dlo_like"], color=colors[mut], s=50, alpha=0.7, zorder=3)
        ax.scatter(
            [fold],
            [sub["frac_1dlo_like"].mean()],
            color=colors[mut],
            s=120,
            marker="D",
            edgecolors="black",
            linewidths=0.8,
            zorder=4,
            label=f"{mut} ({_format_fold(fold)})" if np.isfinite(fold) else mut,
        )
    ax.set_xscale("log")
    criterion = "<" if one_dlo_is_lower else ">"
    ax.set_xlabel("Fold resistance (log scale)")
    ax.set_ylabel(f"Fraction frames 1DLO-like ({criterion} {threshold:.2f} Å)")
    ax.set_title("T290-I63 apo-like tendency vs. DOR resistance", fontweight="bold")
    ax.legend(fontsize=8, frameon=False)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    out_vs_fold = PLOTS / "apo_t290_i63_distance_vs_fold.png"
    fig.savefig(out_vs_fold, dpi=180, bbox_inches="tight")
    plt.close(fig)

def main() -> None:
    ref_1dlo = _reference_distance(REPO / "data/structures/1DLO.cif", I63_CANON, T290_CANON, ATOM_NAME)
    ref_4ncg = _reference_distance(REPO / "data/structures/4NCG.cif", I63_CANON, T290_CANON, ATOM_NAME)
    ref_df = pd.DataFrame(
        [
            {"structure": "1DLO", "leg": "apo", "chain": "p66", "atom_pair": f"{ATOM_NAME}-{ATOM_NAME}", "distance_angstrom": ref_1dlo},
            {"structure": "4NCG", "leg": "holo", "chain": "p66", "atom_pair": f"{ATOM_NAME}-{ATOM_NAME}", "distance_angstrom": ref_4ncg},
        ]
    )
    ref_df.to_csv(OUT_REF, index=False)
    print(f"Reference distances (T290-I63 {ATOM_NAME}-{ATOM_NAME}):")
    for _, r in ref_df.iterrows():
        print(f"  {r['structure']} ({r['leg']}): {r['distance_angstrom']:.3f} Å")
    print(f"Wrote {OUT_REF}")

    all_rows: list[dict] = []
    manifests = [("holo", HOLO_MANIFEST), ("apo", APO_MANIFEST)]
    for leg, manifest in manifests:
        if not manifest.exists():
            continue
        df = pd.read_csv(manifest)
        print(f"Processing {leg} manifest: {manifest} ({len(df)} rows)")
        for _, row in df.iterrows():
            all_rows.extend(_process_replicate(row, leg=leg, frame_stride=1))

    dyn = pd.DataFrame(all_rows)
    if dyn.empty:
        print("No trajectory rows were produced.")
        return
    dyn.to_csv(OUT_DYNAMICS, index=False)
    print(f"Wrote {OUT_DYNAMICS} ({len(dyn)} rows)")

    summary = dyn.groupby(["leg", "mutation", "replicate"], as_index=False).agg(
        mean=("t290_i63_ca_distance_angstrom", "mean"),
        std=("t290_i63_ca_distance_angstrom", "std"),
        count=("t290_i63_ca_distance_angstrom", "count"),
    )
    summary.to_csv(OUT_SUMMARY, index=False)
    print(f"Wrote {OUT_SUMMARY} ({len(summary)} rows)")

    apo = dyn[dyn["leg"] == "apo"].copy()
    apo.to_csv(OUT_APO, index=False)
    print(f"Wrote {OUT_APO} ({len(apo)} rows)")

    fold_map = _load_fold_map()
    _plot_apo(apo, fold_map, ref_1dlo=ref_1dlo, ref_4ncg=ref_4ncg)
    print(f"Wrote {PLOTS / 'apo_t290_i63_distance_timeseries.png'}")
    print(f"Wrote {PLOTS / 'apo_t290_i63_distance_distribution.png'}")
    print(f"Wrote {PLOTS / 'apo_t290_i63_distance_vs_fold.png'}")


if __name__ == "__main__":
    main()
