from __future__ import annotations

import argparse
import json
from pathlib import Path

from .convergence import diagnose_leg, diagnose_phase


def main() -> int:
    parser = argparse.ArgumentParser(description="Diagnose FEP window convergence")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--leg-dir", type=Path)
    group.add_argument("--phase-dir", type=Path)
    parser.add_argument("--target-samples", type=int, default=1000)
    args = parser.parse_args()
    if args.leg_dir:
        result = diagnose_leg(args.leg_dir, target_samples=args.target_samples)
    else:
        result = diagnose_phase(args.phase_dir, target_samples=args.target_samples)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
