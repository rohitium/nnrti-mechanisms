#!/usr/bin/env python3
"""Prepare pmx hybrid protein structures for NEQ mutation FEP legs."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from .alchemical import MutationSite, resolve_mutation_site
from nnrti.fep.mutations import Mutation, MutationLeg, unique_manuscript_legs
from .structure import extract_protein_only, normalize_openmm_for_pmx
from nnrti.fep.config import FEP_PMX_ROOT, LIGAND_RESNAME, PMX_FORCE_FIELD
from nnrti.fep.gromacs_utils import normalize_hybrid_his_for_pdb2gmx


def _leg_by_id(leg_id: str) -> MutationLeg:
    for leg in unique_manuscript_legs():
        if leg.leg_id == leg_id:
            return leg
    raise ValueError(f"Unknown leg_id: {leg_id}")


def _source_pdb(leg: MutationLeg, phase: str, replicate: int) -> Path:
    if phase == "holo":
        return leg.input_complex_pdb(replicate)
    return leg.input_apo_pdb(replicate)


def _endpoint_pdb(leg: MutationLeg, phase: str, replicate: int) -> Path:
    if phase == "holo":
        return leg.endpoint_complex_pdb(replicate)
    return leg.endpoint_apo_pdb(replicate)


def _backend_residue_map(leg_id: str) -> dict | None:
    backend = Path("results/analysis/fep_jorgensen/legs") / leg_id / "prepare_backend.json"
    if not backend.is_file():
        return None
    return json.loads(backend.read_text())


def _resolve_site(leg: MutationLeg, phase: str, replicate: int, chain_id: str) -> MutationSite:
    mutation = Mutation.parse(leg.mutation)
    backend = _backend_residue_map(leg.leg_id)
    if backend and backend.get("pdb_residue_id"):
        return MutationSite(
            chain_id=chain_id,
            ordinal=-1,
            pdb_residue_id=str(backend["pdb_residue_id"]),
            old_residue=mutation.old_residue,
            new_residue=mutation.new_residue,
            mutation=mutation.label,
        )

    start = _source_pdb(leg, phase, replicate)
    end = _endpoint_pdb(leg, phase, replicate)
    if not start.is_file():
        raise FileNotFoundError(start)
    if not end.is_file():
        raise FileNotFoundError(end)
    return resolve_mutation_site(start, end, mutation, chain_id=chain_id)


def _write_mutation_script(path: Path, chain_id: str, pdb_residue_id: str, new_residue: str) -> None:
    path.write_text(f"{chain_id} {pdb_residue_id} {new_residue}\n")


def _run_pmx_mutate(
    *,
    input_pdb: Path,
    output_pdb: Path,
    script_path: Path,
    force_field: str,
    gmxlib: str | None,
    endpoint_pdb: Path | None = None,
) -> None:
    env = os.environ.copy()
    if gmxlib:
        env["GMXLIB"] = gmxlib
    cmd = [
        "pmx",
        "mutate",
        "-f",
        str(input_pdb),
        "-o",
        str(output_pdb),
        "-ff",
        force_field,
        "--script",
        str(script_path),
        "--keep_resid",
    ]
    if endpoint_pdb is not None:
        cmd.extend(["-fB", str(endpoint_pdb)])
    subprocess.run(cmd, check=True, env=env)


def prepare_hybrid(
    leg: MutationLeg,
    *,
    phase: str,
    replicate: int = 1,
    chain_id: str = "A",
    output_dir: Path | None = None,
    force_field: str = PMX_FORCE_FIELD,
    gmxlib: str | None = None,
    run_pdb2gmx: bool = False,
) -> Path:
    """Create pmx hybrid PDB + residue map for one leg/phase/replicate."""
    out_dir = output_dir or (FEP_PMX_ROOT / "legs" / leg.leg_id / phase / f"rep_{replicate:02d}")
    out_dir.mkdir(parents=True, exist_ok=True)

    site = _resolve_site(leg, phase, replicate, chain_id)
    source = _source_pdb(leg, phase, replicate)
    endpoint = _endpoint_pdb(leg, phase, replicate)
    protein_pdb = out_dir / "protein_input.pdb"
    extract_protein_only(
        source,
        out_dir,
        ligand_resname=LIGAND_RESNAME,
        output_name=protein_pdb.name,
    )
    normalize_openmm_for_pmx(protein_pdb)

    endpoint_protein_pdb: Path | None = None
    if endpoint.is_file():
        endpoint_protein_pdb = out_dir / "endpoint_protein.pdb"
        extract_protein_only(
            endpoint,
            out_dir,
            ligand_resname=LIGAND_RESNAME,
            output_name=endpoint_protein_pdb.name,
        )
        normalize_openmm_for_pmx(endpoint_protein_pdb)

    script_path = out_dir / "mutation.script"
    _write_mutation_script(script_path, chain_id, site.pdb_residue_id, site.new_residue)
    hybrid_pdb = out_dir / "hybrid.pdb"
    _run_pmx_mutate(
        input_pdb=protein_pdb,
        output_pdb=hybrid_pdb,
        script_path=script_path,
        force_field=force_field,
        gmxlib=gmxlib,
        endpoint_pdb=endpoint_protein_pdb,
    )
    normalize_hybrid_his_for_pdb2gmx(hybrid_pdb)

    residue_map = {
        "leg_id": leg.leg_id,
        "phase": phase,
        "replicate": replicate,
        "auth_mutation": leg.mutation,
        "chain_id": chain_id,
        "auth_residue_id": Mutation.parse(leg.mutation).residue_id,
        "pdb_residue_id": site.pdb_residue_id,
        "old_residue": site.old_residue,
        "new_residue": site.new_residue,
        "pmx_force_field": force_field,
        "source_complex_pdb": str(source),
        "protein_input_pdb": str(protein_pdb),
        "endpoint_protein_pdb": str(endpoint_protein_pdb) if endpoint_protein_pdb else None,
        "hybrid_pdb": str(hybrid_pdb),
    }
    (out_dir / "residue_map.json").write_text(json.dumps(residue_map, indent=2) + "\n")

    if run_pdb2gmx:
        gmx = shutil.which("gmx") or shutil.which("gmx_mpi")
        if gmx is None:
            raise RuntimeError("gmx not found; load GROMACS before --run-pdb2gmx")
        conf_pdb = out_dir / "conf.pdb"
        topol_top = out_dir / "topol.top"
        subprocess.run(
            [
                gmx,
                "pdb2gmx",
                "-f",
                str(hybrid_pdb),
                "-o",
                str(conf_pdb),
                "-p",
                str(topol_top),
                "-ff",
                force_field,
                "-water",
                "tip3p",
                "-ignh",
            ],
            check=True,
            cwd=out_dir,
        )
        subprocess.run(
            [
                "pmx",
                "gentop",
                "-p",
                str(topol_top),
                "-o",
                str(out_dir / "hybrid.top"),
                "-ff",
                force_field,
            ],
            check=True,
            env={**os.environ, "GMXLIB": gmxlib} if gmxlib else os.environ,
            cwd=out_dir,
        )

    return hybrid_pdb


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Prepare pmx hybrid structures for NEQ FEP.")
    parser.add_argument("--leg", required=True, help="Leg id, e.g. wt_to_V106A")
    parser.add_argument("--phase", choices=("holo", "apo"), default="holo")
    parser.add_argument("--replicate", type=int, default=1)
    parser.add_argument("--chain-id", default="A")
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--gmxlib", default=os.environ.get("GMXLIB", ""))
    parser.add_argument(
        "--run-pdb2gmx",
        action="store_true",
        help="Also run gmx pdb2gmx + pmx gentop (requires GROMACS on PATH).",
    )
    args = parser.parse_args(argv)

    leg = _leg_by_id(args.leg)
    gmxlib = args.gmxlib or None
    hybrid = prepare_hybrid(
        leg,
        phase=args.phase,
        replicate=args.replicate,
        chain_id=args.chain_id,
        output_dir=args.output_dir,
        gmxlib=gmxlib,
        run_pdb2gmx=args.run_pdb2gmx,
    )
    print(f"Wrote hybrid structure: {hybrid}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
