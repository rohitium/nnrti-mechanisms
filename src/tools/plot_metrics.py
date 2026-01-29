from __future__ import annotations

import os
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

# Avoid font cache issues on locked-down home directories.
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".mplconfig"))
os.environ.setdefault("FONTCONFIG_PATH", str(ROOT / ".fontconfig"))

from src.plotting import plot_delta_metrics
from src.utils import ensure_dirs, project_paths


def main() -> None:
    root = Path.cwd()
    paths = project_paths(root)
    ensure_dirs([paths.plots, Path(os.environ["MPLCONFIGDIR"])])
    df = pd.read_csv(paths.results / "metrics_summary.csv")
    plot_delta_metrics(df, paths)
    print(f"Wrote plots to {paths.plots}")


if __name__ == "__main__":
    main()
