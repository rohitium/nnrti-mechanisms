from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import HuberRegressor, LinearRegression, TheilSenRegressor
from sklearn.metrics import r2_score


OUT_DIR = Path("results/analysis/binding_energy/tables")


def _perm_p_corr(x: np.ndarray, y: np.ndarray, stat_fn, rng: np.random.Generator, n: int = 50_000) -> tuple[float, float]:
    obs = float(stat_fn(x, y))
    vals = np.empty(n, dtype=float)
    for i in range(n):
        vals[i] = float(stat_fn(x, rng.permutation(y)))
    p = float((np.sum(np.abs(vals) >= abs(obs)) + 1) / (n + 1))
    return obs, p


def _boot_ci(x: np.ndarray, y: np.ndarray, stat_fn, rng: np.random.Generator, n: int = 20_000) -> tuple[float, float, float]:
    vals: list[float] = []
    m = len(x)
    for _ in range(n):
        idx = rng.integers(0, m, m)
        if len(np.unique(x[idx])) < 2 or len(np.unique(y[idx])) < 2:
            continue
        vals.append(float(stat_fn(x[idx], y[idx])))
    arr = np.asarray(vals, dtype=float)
    lo, mid, hi = np.nanpercentile(arr, [2.5, 50, 97.5])
    return float(lo), float(mid), float(hi)


def _pearson(x: np.ndarray, y: np.ndarray) -> float:
    return float(stats.pearsonr(x, y).statistic)


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    return float(stats.spearmanr(x, y).statistic)


def _kendall(x: np.ndarray, y: np.ndarray) -> float:
    return float(stats.kendalltau(x, y).statistic)


def _model_for_name(name: str):
    if name == "ordinary_least_squares":
        return LinearRegression()
    if name == "huber":
        return HuberRegressor(epsilon=1.35, alpha=0.0)
    if name == "theil_sen":
        return TheilSenRegressor(random_state=1)
    raise ValueError(name)


def main() -> int:
    rng = np.random.default_rng(20260514)
    summary = pd.read_csv(OUT_DIR / "mutation_ddg_summary.csv")
    rep = pd.read_csv(OUT_DIR / "ddg_full.csv")

    rows: list[dict[str, object]] = []
    df = summary.dropna(subset=["fold_reduction", "ddg_electrostatic_mean"]).copy()
    df["x"] = np.log10(pd.to_numeric(df["fold_reduction"], errors="coerce"))
    df["y"] = pd.to_numeric(df["ddg_electrostatic_mean"], errors="coerce")
    df = df.dropna(subset=["x", "y"]).reset_index(drop=True)

    tests = [
        ("pearson", _pearson, stats.pearsonr),
        ("spearman", _spearman, stats.spearmanr),
        ("kendall", _kendall, stats.kendalltau),
    ]
    x = df["x"].to_numpy(dtype=float)
    y = df["y"].to_numpy(dtype=float)
    for name, stat_fn, scipy_fn in tests:
        res = scipy_fn(x, y)
        _obs, pperm = _perm_p_corr(x, y, stat_fn, rng)
        lo, mid, hi = _boot_ci(x, y, stat_fn, rng)
        rows.append(
            {
                "dataset": "mutation_summary_all",
                "method": name,
                "n": len(df),
                "statistic": float(res.statistic),
                "pvalue": float(res.pvalue),
                "permutation_p": pperm,
                "bootstrap_ci95_low": lo,
                "bootstrap_median": mid,
                "bootstrap_ci95_high": hi,
            }
        )

    xm = df[["x"]].to_numpy(dtype=float)
    for name in ["ordinary_least_squares", "huber", "theil_sen"]:
        model = _model_for_name(name)
        model.fit(xm, y)
        pred = model.predict(xm)
        slope = float(model.coef_[0])
        intercept = float(model.intercept_)
        r2 = float(r2_score(y, pred))
        obs = abs(slope)
        vals = np.empty(20_000, dtype=float)
        for i in range(vals.size):
            m = _model_for_name(name)
            m.fit(xm, rng.permutation(y))
            vals[i] = abs(float(m.coef_[0]))
        pperm = float((np.sum(vals >= obs) + 1) / (vals.size + 1))
        rows.append(
            {
                "dataset": "mutation_summary_all",
                "method": name,
                "n": len(df),
                "statistic": slope,
                "pvalue": np.nan,
                "permutation_p": pperm,
                "bootstrap_ci95_low": np.nan,
                "bootstrap_median": intercept,
                "bootstrap_ci95_high": r2,
            }
        )

    loo: list[dict[str, object]] = []
    for _, row in df.iterrows():
        sub = df[df["mutation"] != row["mutation"]]
        loo.append(
            {
                "left_out": row["mutation"],
                "pearson_r": float(stats.pearsonr(sub["x"], sub["y"]).statistic),
                "pearson_p": float(stats.pearsonr(sub["x"], sub["y"]).pvalue),
                "spearman_r": float(stats.spearmanr(sub["x"], sub["y"]).statistic),
                "spearman_p": float(stats.spearmanr(sub["x"], sub["y"]).pvalue),
            }
        )

    sensitivity = [
        ("single_mutations_only", ~df["mutation"].str.contains(r"\+", regex=True)),
        ("combos_only", df["mutation"].str.contains(r"\+", regex=True)),
        ("exclude_y188l", df["mutation"] != "Y188L"),
        ("exclude_v106i_f227c", df["mutation"] != "V106I+F227C"),
        ("exclude_controls_fold_lt2", df["fold_reduction"] >= 2.0),
    ]
    for label, mask in sensitivity:
        sub = df[mask].copy()
        if len(sub) < 5:
            continue
        sx = sub["x"].to_numpy(dtype=float)
        sy = sub["y"].to_numpy(dtype=float)
        for name, stat_fn, scipy_fn in tests[:2]:
            res = scipy_fn(sx, sy)
            _obs, pperm = _perm_p_corr(sx, sy, stat_fn, rng)
            rows.append(
                {
                    "dataset": label,
                    "method": name,
                    "n": len(sub),
                    "statistic": float(res.statistic),
                    "pvalue": float(res.pvalue),
                    "permutation_p": pperm,
                    "bootstrap_ci95_low": np.nan,
                    "bootstrap_median": np.nan,
                    "bootstrap_ci95_high": np.nan,
                }
            )

    rdf = rep.dropna(subset=["fold_reduction", "ddg_electrostatic"]).copy()
    rdf = rdf[rdf["mutation"] != "WT"].copy()
    rdf["x"] = np.log10(pd.to_numeric(rdf["fold_reduction"], errors="coerce"))
    rdf["y"] = pd.to_numeric(rdf["ddg_electrostatic"], errors="coerce")
    rdf = rdf.dropna(subset=["x", "y"])
    rx = rdf["x"].to_numpy(dtype=float)
    ry = rdf["y"].to_numpy(dtype=float)
    for name, _stat_fn, scipy_fn in tests:
        res = scipy_fn(rx, ry)
        rows.append(
            {
                "dataset": "replicate_level_unpooled",
                "method": name,
                "n": len(rdf),
                "statistic": float(res.statistic),
                "pvalue": float(res.pvalue),
                "permutation_p": np.nan,
                "bootstrap_ci95_low": np.nan,
                "bootstrap_median": np.nan,
                "bootstrap_ci95_high": np.nan,
            }
        )

    out = pd.DataFrame(rows)
    loo_df = pd.DataFrame(loo).sort_values("pearson_p")
    out_path = OUT_DIR / "electrostatics_signal_sensitivity_stats.csv"
    loo_path = OUT_DIR / "electrostatics_leave_one_out_stats.csv"
    out.to_csv(out_path, index=False)
    loo_df.to_csv(loo_path, index=False)
    print(out.to_string(index=False))
    print("\nLeave-one-out sorted by Pearson p:")
    print(loo_df.to_string(index=False))
    print(f"\nWrote {out_path}")
    print(f"Wrote {loo_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
