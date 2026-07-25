from __future__ import annotations

import argparse
import json
from pathlib import Path

from .alchemical import build_alchemical_plan, write_holo_phase
from .config import FEPConfig
from .mutations import MutationLeg
from .perses_hybrid import (
    perses_available,
    prepare_apo_hybrid,
    prepare_hybrid_leg,
    prepare_holo_hybrid,
)

PHASES = ("holo", "apo", "all")


def _normalize_phases(value: str) -> tuple[str, ...]:
    if value == "all":
        return ("holo", "apo")
    if value in {"holo", "apo"}:
        return (value,)
    raise ValueError(f"Unknown phase selection: {value!r}")


def prepare(
    config: FEPConfig,
    replicate: int = 1,
    phases: tuple[str, ...] = ("holo", "apo"),
) -> None:
    """Prepare holo and/or apo mutation legs for MBAR sampling."""
    backend = config.prepare_backend
    selected = tuple(dict.fromkeys(phases))
    if backend == "perses":
        if not perses_available():
            raise ImportError(
                "Perses hybrid prep requires perses, openmmtools, ambertools GAFF, RDKit, and "
                "OpenFF toolkit. Run: bash scripts/fep_jorgensen/setup_perses_env.sh "
                "or use --backend scaling."
            )
        if selected == ("holo",):
            prepare_holo_hybrid(config, replicate=replicate)
        elif selected == ("apo",):
            prepare_apo_hybrid(config, replicate=replicate)
        else:
            prepare_hybrid_leg(config, replicate=replicate, phases=selected)
        return
    if backend != "scaling":
        raise ValueError(f"Unknown prepare_backend: {backend!r}")
    if selected != ("holo",):
        raise ValueError("Scaling backend currently supports holo phase only")

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
    (config.run_dir / "prepare_backend.json").write_text(
        json.dumps({"backend": "scaling", "phases": ["holo"]}, indent=2) + "\n"
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Prepare holo and/or apo mutation FEP legs"
    )
    parser.add_argument("--mutation", default="V106A", help="Single substitution made in this leg")
    parser.add_argument("--start-label", default="WT")
    parser.add_argument("--end-label", help="Resulting single or compound mutant label")
    parser.add_argument("--input-complex-pdb", "--wt-complex-pdb", dest="input_complex_pdb", type=Path)
    parser.add_argument("--replicate", type=int, default=1)
    parser.add_argument("--output-dir", type=Path, default=FEPConfig().output_dir)
    parser.add_argument(
        "--phase",
        choices=PHASES,
        default="all",
        help="Thermodynamic phase(s) to prepare (default: holo + apo)",
    )
    parser.add_argument(
        "--backend",
        choices=("perses", "scaling"),
        default=FEPConfig().prepare_backend,
        help="perses: hybrid topology; scaling: MD-asset nonbonded scaling fallback",
    )
    args = parser.parse_args()
    end_label = args.end_label or args.mutation
    leg = MutationLeg(args.start_label, end_label, args.mutation)
    overrides = {"output_dir": args.output_dir, "prepare_backend": args.backend}
    if args.input_complex_pdb:
        overrides["wt_complex_pdb"] = args.input_complex_pdb
    config = FEPConfig.for_leg(leg, **overrides)
    prepare(config, replicate=args.replicate, phases=_normalize_phases(args.phase))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
