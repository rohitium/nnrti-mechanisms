#!/usr/bin/env python3
"""V106A worked-example FEP protocol figures (thin wrapper).

Prefer the multi-genotype entry point:

    python3 scripts/fep_pmx/plot_protocol_figures.py --targets V106A

This wrapper regenerates only V106A (also mirrors to protocol_v106a/).
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.fep_pmx.plot_protocol_figures import main


if __name__ == "__main__":
    raise SystemExit(main(["--targets", "V106A", *sys.argv[1:]]))
