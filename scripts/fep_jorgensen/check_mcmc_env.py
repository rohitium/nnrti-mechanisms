"""Verify imports required for Perses/openmmtools MCMC sampling."""

from __future__ import annotations

import json
import sys


def main() -> int:
    report: dict[str, str | bool] = {"ok": True}
    try:
        from scripts.fep_jorgensen.perses_hybrid import perses_available

        report["perses_available"] = perses_available()
        if not report["perses_available"]:
            report["ok"] = False
            report["error"] = "perses_available() returned False"
            print(json.dumps(report, indent=2))
            return 1

        import numpy  # noqa: F401
        import openmm  # noqa: F401
        import openmmtools  # noqa: F401
        import perses  # noqa: F401
        from openmmtools import mcmc  # noqa: F401
        from openmmtools.multistate import MultiStateReporter  # noqa: F401
        from perses.samplers.multistate import HybridRepexSampler  # noqa: F401

        report["numpy"] = numpy.__version__
        report["openmm"] = openmm.__version__
        report["openmmtools"] = openmmtools.__version__
        report["perses"] = perses.__version__
    except Exception as exc:
        report["ok"] = False
        report["error"] = str(exc)

    print(json.dumps(report, indent=2))
    return 0 if report["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
