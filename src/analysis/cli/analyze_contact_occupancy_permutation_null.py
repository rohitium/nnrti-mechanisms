from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a trajectory-label permutation null for WT-referenced DOR contact "
            "occupancy shifts."
        )
    )
    parser.add_argument(
        "--replicate-contact-csv",
        type=Path,
        default=Path("results/analysis/triplet_story_analyses/contact_story_all_mutations_excluding_f227c/tables/replicate_contact.csv"),
    )
    parser.add_argument(
        "--display-residue-csv",
        type=Path,
        default=Path(
            "results/analysis/triplet_story_analyses/contact_story_all_mutations_excluding_f227c/tables/"
            "all_mutation_wt_referenced_occupancy_heatmap_wt_contacted_residues_by_region_excluding_f227c_display_residues.csv"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/analysis/triplet_story_analyses/contact_occupancy_permutation_null"),
    )
    parser.add_argument("--n-permutations", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=20260517)
    parser.add_argument("--alpha", type=float, default=0.10)
    return parser.parse_args()


def _complete_contact_grid(rep: pd.DataFrame, residue_df: pd.DataFrame) -> tuple[pd.DataFrame, list[tuple[int, int]]]:
    keys = list(zip(residue_df["traj_resid"].astype(int), residue_df["auth_resid"].astype(int)))
    rep = rep.copy()
    rep["traj_resid"] = pd.to_numeric(rep["traj_resid"], errors="coerce").astype("Int64")
    rep["auth_resid"] = pd.to_numeric(rep["auth_resid"], errors="coerce").astype("Int64")
    rep["replicate"] = pd.to_numeric(rep["replicate"], errors="coerce").astype("Int64")
    rep["key"] = list(zip(rep["traj_resid"].astype(int), rep["auth_resid"].astype(int)))
    rep["traj_id"] = rep["mutation"].astype(str) + "_rep" + rep["replicate"].astype(str)

    frames = (
        rep.groupby(["mutation", "replicate", "traj_id"], as_index=False)["n_total_frames"]
        .max()
        .sort_values(["mutation", "replicate"], kind="stable")
    )
    grid = pd.MultiIndex.from_product([frames["traj_id"].tolist(), keys], names=["traj_id", "key"]).to_frame(index=False)
    grid = grid.merge(frames[["traj_id", "mutation", "replicate", "n_total_frames"]], on="traj_id", how="left")
    observed = rep[rep["key"].isin(keys)][["traj_id", "key", "n_contact_frames"]].copy()
    out = grid.merge(observed, on=["traj_id", "key"], how="left")
    out["n_contact_frames"] = out["n_contact_frames"].fillna(0).astype(int)
    return out, keys


def _matrix_from_grid(grid: pd.DataFrame, keys: list[tuple[int, int]]) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    traj_df = (
        grid[["traj_id", "mutation", "replicate", "n_total_frames"]]
        .drop_duplicates()
        .sort_values(["mutation", "replicate"], kind="stable")
        .reset_index(drop=True)
    )
    traj_index = {traj_id: i for i, traj_id in enumerate(traj_df["traj_id"].tolist())}
    key_index = {key: i for i, key in enumerate(keys)}

    contacts = np.zeros((len(traj_df), len(keys)), dtype=float)
    totals = np.repeat(traj_df["n_total_frames"].to_numpy(dtype=float)[:, None], len(keys), axis=1)
    for row in grid.itertuples(index=False):
        contacts[traj_index[row.traj_id], key_index[row.key]] = float(row.n_contact_frames)
    return contacts, totals, traj_df["mutation"].to_numpy(dtype=str), traj_df["traj_id"].tolist()


def _pooled_means(
    contacts: np.ndarray,
    totals: np.ndarray,
    labels: np.ndarray,
    mutation_order: list[str],
) -> dict[str, np.ndarray]:
    out: dict[str, np.ndarray] = {}
    for mutation in mutation_order:
        mask = labels == mutation
        if not np.any(mask):
            raise ValueError(f"No trajectories found for mutation label {mutation!r}")
        out[mutation] = contacts[mask].sum(axis=0) / totals[mask].sum(axis=0)
    return out


def _plot_summary(summary: pd.DataFrame, output_png: Path, global_threshold: float, threshold_percentile: float) -> None:
    import matplotlib.pyplot as plt

    plot_df = summary.sort_values("max_abs_shift", ascending=True).copy()
    colors = np.where(plot_df["exceeds_global_threshold"], "#d04f45", "#5f7f9f")
    fig_h = max(6.5, 0.34 * len(plot_df) + 1.6)
    fig, ax = plt.subplots(figsize=(10.5, fig_h), constrained_layout=True)
    y = np.arange(len(plot_df))
    ax.barh(y, plot_df["max_abs_shift"], color=colors, edgecolor="#333333", linewidth=0.5)
    ax.axvline(
        global_threshold,
        color="#333333",
        linestyle="--",
        linewidth=1.5,
    )
    ax.set_yticks(y, plot_df["mutation"].tolist())
    ax.set_xlabel("|WT-referenced occupancy shift|", fontsize=18, fontweight="bold")
    ax.set_ylabel("")
    ax.tick_params(axis="both", labelsize=14)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", alpha=0.24, linestyle=":")
    ax.text(
        global_threshold + 0.012,
        -0.62,
        f"Global null {threshold_percentile:.0f}th percentile",
        ha="left",
        va="center",
        fontsize=13,
        color="#333333",
    )
    for yi, row in enumerate(plot_df.itertuples(index=False)):
        ax.text(
            float(row.max_abs_shift) + 0.012,
            yi,
            f"{row.top_residue}",
            va="center",
            ha="left",
            fontsize=11,
        )
    ax.set_xlim(0, max(float(plot_df["max_abs_shift"].max()), global_threshold) * 1.18)
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=320, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = _parse_args()
    if int(args.n_permutations) < 100:
        raise ValueError("--n-permutations should be at least 100 for a stable empirical null.")
    if not (0 < float(args.alpha) < 1):
        raise ValueError("--alpha must be between 0 and 1.")

    out_tables = args.output_dir / "tables"
    out_plots = args.output_dir / "plots"
    out_config = args.output_dir / "config"
    for d in (out_tables, out_plots, out_config):
        d.mkdir(parents=True, exist_ok=True)

    rep = pd.read_csv(args.replicate_contact_csv)
    residue_df = pd.read_csv(args.display_residue_csv)
    grid, keys = _complete_contact_grid(rep, residue_df)
    input_keyed = rep.copy()
    input_keyed["key"] = list(zip(input_keyed["traj_resid"].astype(int), input_keyed["auth_resid"].astype(int)))
    n_observed_grid_rows = int(len(input_keyed[input_keyed["key"].isin(keys)]))
    contacts, totals, true_labels, traj_ids = _matrix_from_grid(grid, keys)
    mutation_order = ["WT"] + [m for m in sorted(pd.unique(true_labels)) if m != "WT"]
    mutant_order = [m for m in mutation_order if m != "WT"]

    label_by_key = dict(zip(keys, residue_df["label"].astype(str)))
    region_by_key = dict(zip(keys, residue_df.get("pocket_region", pd.Series([""] * len(residue_df))).astype(str)))

    observed_means = _pooled_means(contacts, totals, true_labels, mutation_order)
    observed_rows: list[dict[str, object]] = []
    for mutation in mutant_order:
        shift = observed_means[mutation] - observed_means["WT"]
        for key, value in zip(keys, shift):
            observed_rows.append(
                {
                    "mutation": mutation,
                    "residue": label_by_key[key],
                    "pocket_region": region_by_key.get(key, ""),
                    "traj_resid": int(key[0]),
                    "auth_resid": int(key[1]),
                    "wt_referenced_occupancy_shift": float(value),
                    "abs_wt_referenced_occupancy_shift": float(abs(value)),
                }
            )
    observed_df = pd.DataFrame(observed_rows)

    rng = np.random.default_rng(int(args.seed))
    n_perm = int(args.n_permutations)
    null_global_max = np.empty(n_perm, dtype=float)
    null_mutation_max = {mutation: np.empty(n_perm, dtype=float) for mutation in mutant_order}

    for i in range(n_perm):
        permuted_labels = rng.permutation(true_labels)
        permuted_means = _pooled_means(contacts, totals, permuted_labels, mutation_order)
        mutation_maxima: list[float] = []
        for mutation in mutant_order:
            max_shift = float(np.max(np.abs(permuted_means[mutation] - permuted_means["WT"])))
            null_mutation_max[mutation][i] = max_shift
            mutation_maxima.append(max_shift)
        null_global_max[i] = max(mutation_maxima)

    threshold_percentile = 100.0 * (1.0 - float(args.alpha))
    global_threshold = float(np.quantile(null_global_max, 1.0 - float(args.alpha)))

    top_rows: list[dict[str, object]] = []
    for mutation, g in observed_df.groupby("mutation", sort=False):
        top = g.sort_values("abs_wt_referenced_occupancy_shift", ascending=False).iloc[0]
        max_abs = float(top["abs_wt_referenced_occupancy_shift"])
        mutation_null = null_mutation_max[str(mutation)]
        top_rows.append(
            {
                "mutation": mutation,
                "top_residue": str(top["residue"]),
                "top_residue_region": str(top["pocket_region"]),
                "top_wt_referenced_occupancy_shift": float(top["wt_referenced_occupancy_shift"]),
                "max_abs_shift": max_abs,
                "mutation_null_threshold_percentile": float(threshold_percentile),
                "mutation_null_threshold": float(np.quantile(mutation_null, 1.0 - float(args.alpha))),
                "global_null_threshold_percentile": float(threshold_percentile),
                "global_null_threshold": global_threshold,
                "mutation_level_empirical_p": float((np.sum(mutation_null >= max_abs) + 1) / (n_perm + 1)),
                "global_fwer_p": float((np.sum(null_global_max >= max_abs) + 1) / (n_perm + 1)),
                "exceeds_global_threshold": bool(max_abs > global_threshold),
            }
        )
    summary_df = pd.DataFrame(top_rows).sort_values("max_abs_shift", ascending=False)

    observed_df["global_fwer_p_for_cell"] = observed_df["abs_wt_referenced_occupancy_shift"].map(
        lambda x: float((np.sum(null_global_max >= float(x)) + 1) / (n_perm + 1))
    )
    observed_df["exceeds_global_threshold"] = observed_df["abs_wt_referenced_occupancy_shift"] > global_threshold
    candidates_df = observed_df[observed_df["exceeds_global_threshold"]].sort_values(
        "abs_wt_referenced_occupancy_shift",
        ascending=False,
    )

    null_df = pd.DataFrame({"permutation": np.arange(1, n_perm + 1), "global_max_abs_shift": null_global_max})
    for mutation in mutant_order:
        null_df[f"{mutation}_max_abs_shift"] = null_mutation_max[mutation]

    summary_df.to_csv(out_tables / "mutation_max_shift_permutation_summary.csv", index=False)
    observed_df.sort_values("abs_wt_referenced_occupancy_shift", ascending=False).to_csv(
        out_tables / "residue_shift_permutation_results.csv",
        index=False,
    )
    candidates_df.to_csv(out_tables / "candidate_reporter_coordinates.csv", index=False)
    null_df.to_csv(out_tables / "permutation_null_distribution.csv", index=False)

    _plot_summary(
        summary=summary_df,
        output_png=out_plots / "mutation_max_shift_vs_permutation_null.png",
        global_threshold=global_threshold,
        threshold_percentile=threshold_percentile,
    )

    config = {
        "replicate_contact_csv": str(args.replicate_contact_csv),
        "display_residue_csv": str(args.display_residue_csv),
        "output_dir": str(args.output_dir),
        "n_permutations": n_perm,
        "seed": int(args.seed),
        "alpha": float(args.alpha),
        "n_trajectories": int(len(traj_ids)),
        "n_mutations_including_wt": int(len(mutation_order)),
        "n_mutant_comparisons": int(len(mutant_order)),
        "n_residues_screened": int(len(keys)),
        "zero_contact_rows_added": int(len(grid) - n_observed_grid_rows),
        "global_null_threshold_percentile": threshold_percentile,
        "global_null_threshold": global_threshold,
    }
    (out_config / "run_config.json").write_text(json.dumps(config, indent=2))
    print(f"Saved {out_tables / 'mutation_max_shift_permutation_summary.csv'}")
    print(f"Saved {out_tables / 'candidate_reporter_coordinates.csv'}")
    print(f"Saved {out_plots / 'mutation_max_shift_vs_permutation_null.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
