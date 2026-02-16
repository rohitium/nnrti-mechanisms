#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from .plot_all_mutation_drm_distances import curate_interesting_traces


def _norm_mutation_key(m: str) -> str:
    return str(m).strip().replace("+", "_").upper()


def _enrich_trace_paths(df: pd.DataFrame, manifest_csv: Path) -> pd.DataFrame:
    if not manifest_csv.exists() or df.empty:
        return df
    mf = pd.read_csv(manifest_csv)
    need_cols = {"mutation", "replicate", "output_json", "safe_label"}
    if not need_cols.issubset(set(mf.columns)):
        return df

    map_rows = mf[list(need_cols)].copy()
    map_rows["mutation_key"] = map_rows["mutation"].astype(str).map(_norm_mutation_key)
    map_rows["replicate"] = pd.to_numeric(map_rows["replicate"], errors="coerce").astype("Int64")
    map_rows = map_rows.dropna(subset=["replicate"])
    map_rows["replicate"] = map_rows["replicate"].astype(int)

    json_cache: dict[str, dict] = {}
    def _json_value(path_text: str, key: str) -> str:
        p = str(path_text or "").strip()
        if not p:
            return ""
        if p not in json_cache:
            try:
                json_cache[p] = json.loads(Path(p).read_text())
            except Exception:
                json_cache[p] = {}
        return str(json_cache[p].get(key, "") or "")

    merged = df.copy()
    merged["mutation_key"] = merged["mutation"].astype(str).map(_norm_mutation_key)
    merged["replicate"] = pd.to_numeric(merged["replicate"], errors="coerce").astype("Int64")
    merged = merged.dropna(subset=["replicate"])
    merged["replicate"] = merged["replicate"].astype(int)
    merged = merged.merge(
        map_rows[["mutation_key", "replicate", "output_json", "safe_label"]],
        on=["mutation_key", "replicate"],
        how="left",
        suffixes=("", "_mf"),
    )

    if "output_json_mf" in merged.columns:
        merged["output_json"] = merged["output_json"].fillna("")
        merged.loc[merged["output_json"].astype(str).str.strip() == "", "output_json"] = merged["output_json_mf"]
    if "safe_label_mf" in merged.columns:
        merged["safe_label"] = merged["safe_label"].fillna("")
        merged.loc[merged["safe_label"].astype(str).str.strip() == "", "safe_label"] = merged["safe_label_mf"]

    if "analysis_dcd" in merged.columns:
        merged["analysis_dcd"] = merged["analysis_dcd"].fillna("")
        mask = merged["analysis_dcd"].astype(str).str.strip() == ""
        merged.loc[mask, "analysis_dcd"] = merged.loc[mask, "output_json"].map(
            lambda p: _json_value(str(p), "analysis_dcd")
        )
    if "analysis_topology_pdb" in merged.columns:
        merged["analysis_topology_pdb"] = merged["analysis_topology_pdb"].fillna("")
        mask = merged["analysis_topology_pdb"].astype(str).str.strip() == ""
        merged.loc[mask, "analysis_topology_pdb"] = merged.loc[mask, "output_json"].map(
            lambda p: _json_value(str(p), "analysis_topology_pdb")
        )
    return merged.drop(columns=[c for c in ["mutation_key", "output_json_mf", "safe_label_mf"] if c in merged.columns])


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Curate a ranked list of interesting DRM distance traces from existing timeseries CSV."
    )
    parser.add_argument(
        "--input-csv",
        type=Path,
        default=Path("results/drm_sidechain_distance_timeseries_all_mutations.csv"),
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=Path("results/drm_sidechain_distance_interesting_traces.csv"),
    )
    parser.add_argument(
        "--plots-dir",
        type=Path,
        default=Path("results/plots/drm_distances"),
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("results/md_manifest.csv"),
        help="Optional manifest used to backfill output_json and analysis file paths.",
    )
    parser.add_argument("--top-n", type=int, default=60)
    parser.add_argument("--min-score", type=float, default=5.0)
    args = parser.parse_args()

    if not args.input_csv.exists():
        raise FileNotFoundError(f"Missing input timeseries CSV: {args.input_csv}")

    ts_df = pd.read_csv(args.input_csv)
    out_df = curate_interesting_traces(
        ts_df,
        output_csv=args.output_csv,
        plots_dir=args.plots_dir,
        top_n=int(args.top_n),
        min_score=float(args.min_score),
    )
    out_df = _enrich_trace_paths(out_df, args.manifest)
    out_df.to_csv(args.output_csv, index=False)
    print(f"Wrote {args.output_csv} (rows={len(out_df)}, interesting={int(out_df['is_interesting'].sum()) if not out_df.empty else 0})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
