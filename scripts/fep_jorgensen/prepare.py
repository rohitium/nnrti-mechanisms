from __future__ import annotations

import argparse
import json
from pathlib import Path

from .alchemical import build_alchemical_plan, write_holo_phase
from .config import FEPConfig
from .mutations import MutationLeg


def prepare(config: FEPConfig, replicate: int = 1) -> None:
    """Prepare one holo mutation leg from existing MD assets.

    Uses the start genotype's serialized OpenMM system plus endpoint PDB diffing
    to identify the alchemical side-chain atoms.  No Perses/OpenEye rebuild is
    required because the manuscript MD pipeline already produced equilibrated
    holo complexes for every genotype.
    """
    config.validate(require_inputs=True)
    plan = build_alchemical_plan(config.leg, replicate=replicate, chain_id=config.chain_id)
    write_holo_phase(
        plan,
        config.run_dir / "holo",
        config.lambda_schedule.values,
        config.leg,
    )
    config.write(config.run_dir / "config.json")
    config.approx_protocol.write(config.run_dir / "approx_protocol.json")
    (config.run_dir / "alchemical_plan.json").write_text(
        json.dumps(plan.to_dict(), indent=2) + "\n"
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Prepare one holo mutation FEP leg from existing MD assets"
    )
    parser.add_argument("--mutation", default="V106A", help="Single substitution made in this leg")
    parser.add_argument("--start-label", default="WT")
    parser.add_argument("--end-label", help="Resulting single or compound mutant label")
    parser.add_argument("--input-complex-pdb", "--wt-complex-pdb", dest="input_complex_pdb", type=Path)
    parser.add_argument("--replicate", type=int, default=1)
    parser.add_argument("--output-dir", type=Path, default=FEPConfig().output_dir)
    args = parser.parse_args()
    end_label = args.end_label or args.mutation
    leg = MutationLeg(args.start_label, end_label, args.mutation)
    overrides = {"output_dir": args.output_dir}
    if args.input_complex_pdb:
        overrides["wt_complex_pdb"] = args.input_complex_pdb
    config = FEPConfig.for_leg(leg, **overrides)
    prepare(config, replicate=args.replicate)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
