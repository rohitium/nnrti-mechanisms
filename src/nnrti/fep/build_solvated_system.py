#!/usr/bin/env python3
"""Build solvated GROMACS hybrid systems for pmx NEQ FEP."""

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

from nnrti.fep.config import (
    BOX_TYPE,
    CHARGE_LEG_DELTA_Q,
    DOR_ITP_DIR,
    FEP_PMX_ROOT,
    IONIC_STRENGTH_M,
    LIGAND_RESNAME,
    PMX_FORCE_FIELD_LABEL,
    SOLVENT_PADDING_NM,
    USE_COALCHEMICAL_ION,
    WATER_MODEL,
)
from nnrti.fep.coalchemical_ion import add_coalchemical_ion
from nnrti.fep.gromacs_utils import (
    GromacsError,
    append_ligand_to_topology,
    count_net_charge_from_top,
    extract_ligand_coords_nm,
    find_gmx,
    normalize_hybrid_his_for_pdb2gmx,
    parse_dor_atom_names,
    parse_gro_atom_count,
    run_gmx,
    write_dor_ligand_itps,
    write_ligand_gro,
    write_merged_gro,
)
from nnrti.fep.prepare_hybrid import _leg_by_id, _source_pdb


def _default_hybrid_dir(leg_id: str, phase: str, replicate: int) -> Path:
    return FEP_PMX_ROOT / "legs" / leg_id / phase / f"rep_{replicate:02d}"


def _default_build_dir(leg_id: str, phase: str, replicate: int) -> Path:
    return _default_hybrid_dir(leg_id, phase, replicate) / "gromacs_build"


def build_solvated_system(
    leg_id: str,
    *,
    phase: str,
    replicate: int = 1,
    hybrid_pdb: Path | None = None,
    build_dir: Path | None = None,
    gmxlib: str | None = None,
    dor_dir: Path = DOR_ITP_DIR,
    force_field: str = PMX_FORCE_FIELD_LABEL,
    validate_grompp: bool = True,
    force: bool = False,
) -> Path:
    """Run pdb2gmx → gentop → (holo ligand) → box → solvate → genion."""
    leg = _leg_by_id(leg_id)
    hybrid_dir = _default_hybrid_dir(leg_id, phase, replicate)
    hybrid_pdb = hybrid_pdb or (hybrid_dir / "hybrid.pdb")
    if not hybrid_pdb.is_file():
        raise FileNotFoundError(f"Missing hybrid PDB: {hybrid_pdb}")

    out_dir = build_dir or _default_build_dir(leg_id, phase, replicate)
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = out_dir / "build_manifest.json"

    final_gro = out_dir / "system.gro"
    final_top = out_dir / "system.top"
    if final_gro.is_file() and final_top.is_file() and not force:
        return final_gro

    gmx = find_gmx()
    env = os.environ.copy()
    if gmxlib:
        env["GMXLIB"] = gmxlib

    em_mdp = out_dir / "em.mdp"
    shutil.copy2(REPO_ROOT / "ops/slurm/fep/mdp/em.mdp", em_mdp)

    work_hybrid = out_dir / "hybrid.pdb"
    shutil.copy2(hybrid_pdb, work_hybrid)
    normalize_hybrid_his_for_pdb2gmx(work_hybrid)

    conf_gro = out_dir / "conf.gro"
    topol_top = out_dir / "topol.top"
    hybrid_top = out_dir / "topol_hybrid.top"

    # RT dimer: 4 termini prompts (NH3+/COO- defaults = option 1).
    ter_input = "1\n" * 4
    run_gmx(
        gmx,
        [
            "pdb2gmx",
            "-f",
            "hybrid.pdb",
            "-o",
            "conf.gro",
            "-p",
            "topol.top",
            "-ff",
            force_field,
            "-water",
            WATER_MODEL,
            "-merge",
            "no",
        ],
        cwd=out_dir,
        input_text=ter_input,
        env=env,
    )

    gentop = subprocess.run(
        [
            "pmx",
            "gentop",
            "-p",
            "topol.top",
            "-o",
            "topol_hybrid.top",
            "-ff",
            force_field,
        ],
        cwd=out_dir,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    if gentop.returncode != 0:
        raise GromacsError(f"pmx gentop failed:\n{gentop.stdout}\n{gentop.stderr}")

    solute_gro = conf_gro
    system_top = hybrid_top

    if phase == "holo":
        dor_top = dor_dir / "dor.top"
        if not dor_top.is_file():
            raise FileNotFoundError(f"Missing DOR topology: {dor_top}")

        ligand_atomtypes_itp = out_dir / "dor_atomtypes.itp"
        ligand_itp = out_dir / "dor.itp"
        write_dor_ligand_itps(
            dor_top,
            atomtypes_itp=ligand_atomtypes_itp,
            molecule_itp=ligand_itp,
        )
        atom_names = parse_dor_atom_names(dor_top, LIGAND_RESNAME)
        source_complex = _source_pdb(leg, phase, replicate)
        coords = extract_ligand_coords_nm(
            source_complex,
            resname=LIGAND_RESNAME,
            atom_names=atom_names,
        )
        ligand_gro = out_dir / "ligand.gro"
        write_ligand_gro(
            coords_nm=coords,
            atom_names=atom_names,
            resname=LIGAND_RESNAME,
            output_gro=ligand_gro,
        )
        complex_gro = out_dir / "complex.gro"
        write_merged_gro(
            protein_gro=conf_gro,
            ligand_gro=ligand_gro,
            output_gro=complex_gro,
            title="Hybrid protein + DOR",
        )
        solute_gro = complex_gro

        holo_top = out_dir / "topol_holo.top"
        shutil.copy2(hybrid_top, holo_top)
        append_ligand_to_topology(
            holo_top,
            ligand_atomtypes_itp=ligand_atomtypes_itp,
            ligand_itp=ligand_itp,
            ligand_name=LIGAND_RESNAME,
            output_top=out_dir / "topol_holo_merged.top",
        )
        system_top = out_dir / "topol_holo_merged.top"

    box_gro = out_dir / "box.gro"
    run_gmx(
        gmx,
        [
            "editconf",
            "-f",
            solute_gro.name,
            "-o",
            box_gro.name,
            "-c",
            "-d",
            str(SOLVENT_PADDING_NM),
            "-bt",
            BOX_TYPE,
        ],
        cwd=out_dir,
        env=env,
    )

    solv_gro = out_dir / "solv.gro"
    solv_top = out_dir / "solv.top"
    shutil.copy2(system_top, solv_top)
    run_gmx(
        gmx,
        [
            "solvate",
            "-cp",
            box_gro.name,
            "-cs",
            "spc216.gro",
            "-o",
            solv_gro.name,
            "-p",
            solv_top.name,
        ],
        cwd=out_dir,
        env=env,
    )

    ions_tpr = out_dir / "ions.tpr"
    run_gmx(
        gmx,
        [
            "grompp",
            "-f",
            em_mdp.name,
            "-c",
            solv_gro.name,
            "-p",
            solv_top.name,
            "-o",
            ions_tpr.name,
            "-maxwarn",
            "10",
        ],
        cwd=out_dir,
        env=env,
    )
    run_gmx(
        gmx,
        [
            "genion",
            "-s",
            ions_tpr.name,
            "-o",
            final_gro.name,
            "-p",
            solv_top.name,
            "-pname",
            "NA",
            "-nname",
            "CL",
            "-neutral",
            "-conc",
            str(IONIC_STRENGTH_M),
        ],
        cwd=out_dir,
        input_text="SOL\n",
        env=env,
    )
    shutil.copy2(solv_top, final_top)

    # Charge-changing legs: the co-alchemical ion is ABANDONED (does not converge;
    # see config.USE_COALCHEMICAL_ION). Default is a RAW non-neutral box + analytical
    # net-charge correction post-hoc. The block below only runs if explicitly
    # re-enabled to reproduce the abandoned experiment.
    coalch_info = None
    delta_q = CHARGE_LEG_DELTA_Q.get(leg_id)
    if delta_q is not None and USE_COALCHEMICAL_ION:
        # COALCH_ION_RANK (default 0 = farthest ion) lets the placement-insensitivity
        # check rebuild a leg using a different bulk ion without editing code.
        ion_rank = int(os.environ.get("COALCH_ION_RANK", "0"))
        coalch_info = add_coalchemical_ion(final_top, final_gro, delta_q=delta_q, ion_rank=ion_rank)

    n_atoms = parse_gro_atom_count(final_gro)
    net_charge = count_net_charge_from_top(final_top)

    if validate_grompp:
        run_gmx(
            gmx,
            [
                "grompp",
                "-f",
                em_mdp.name,
                "-c",
                final_gro.name,
                "-p",
                final_top.name,
                "-o",
                "em.tpr",
                "-maxwarn",
                "10",
            ],
            cwd=out_dir,
            env=env,
        )

    manifest = {
        "leg_id": leg_id,
        "phase": phase,
        "replicate": replicate,
        "force_field": force_field,
        "box_type": BOX_TYPE,
        "solvent_padding_nm": SOLVENT_PADDING_NM,
        "ionic_strength_m": IONIC_STRENGTH_M,
        "hybrid_pdb": str(hybrid_pdb),
        "system_gro": str(final_gro),
        "system_top": str(final_top),
        "n_atoms": n_atoms,
        "net_charge_before_ions": net_charge,
        "validated_grompp": validate_grompp,
        "coalchemical_ion": coalch_info,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    return final_gro


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build solvated GROMACS hybrid system for NEQ FEP.")
    parser.add_argument("--leg", required=True, help="Leg id, e.g. wt_to_V106A")
    parser.add_argument("--phase", choices=("holo", "apo"), required=True)
    parser.add_argument("--replicate", type=int, default=1)
    parser.add_argument("--hybrid-pdb", type=Path, default=None)
    parser.add_argument("--build-dir", type=Path, default=None)
    parser.add_argument("--gmxlib", default=os.environ.get("GMXLIB", ""))
    parser.add_argument("--dor-dir", type=Path, default=DOR_ITP_DIR)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--skip-grompp", action="store_true")
    args = parser.parse_args(argv)

    gmxlib = args.gmxlib or None
    try:
        gro = build_solvated_system(
            args.leg,
            phase=args.phase,
            replicate=args.replicate,
            hybrid_pdb=args.hybrid_pdb,
            build_dir=args.build_dir,
            gmxlib=gmxlib,
            dor_dir=args.dor_dir,
            force=args.force,
            validate_grompp=not args.skip_grompp,
        )
    except (GromacsError, FileNotFoundError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"Built solvated system: {gro}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
