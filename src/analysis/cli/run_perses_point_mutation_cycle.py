from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import time


AA1_TO_AA3 = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "Q": "GLN",
    "E": "GLU",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
}

SOLVENT_RESNAMES = {"HOH", "WAT", "SOL"}
ION_RESNAMES = {"NA", "CL", "K", "Na+", "Cl-"}


def parse_mutation_label(label: str) -> tuple[str, str, str]:
    match = re.fullmatch(r"([A-Z])(\d+)([A-Z])", label.strip().upper())
    if not match:
        raise ValueError(f"Expected single-point mutation like V106A, got {label!r}")
    old_aa, resid, new_aa = match.groups()
    return AA1_TO_AA3[old_aa], resid, AA1_TO_AA3[new_aa]


def extract_protein_and_ligand(
    complex_pdb: Path,
    ligand_template_sdf: Path,
    output_dir: Path,
    ligand_resname: str = "2KW",
) -> tuple[Path, Path]:
    """Create Perses inputs from an existing bound RT-DOR PDB.

    The protein input keeps protein/DNA/cofactor atoms and removes water, ions,
    and DOR. The ligand SDF keeps the existing bond order from the template SDF
    but replaces coordinates with the bound DOR coordinates from the PDB.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    protein_pdb = output_dir / "wt_receptor_no_ligand.pdb"
    ligand_sdf = output_dir / "dor_bound_pose.sdf"

    ligand_xyz = []
    protein_lines = []
    for line in complex_pdb.read_text().splitlines():
        record = line[:6].strip()
        if record in {"ATOM", "HETATM"}:
            resname = line[17:20].strip()
            if resname == ligand_resname:
                ligand_xyz.append(
                    (
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54]),
                    )
                )
                continue
            if resname in SOLVENT_RESNAMES or resname in ION_RESNAMES:
                continue
            protein_lines.append(line)
        elif line.startswith(("CRYST1", "MODEL", "ENDMDL")):
            continue
    if not ligand_xyz:
        raise ValueError(f"No {ligand_resname} ligand atoms found in {complex_pdb}")

    protein_pdb.write_text("\n".join(protein_lines + ["END", ""]))

    sdf_lines = ligand_template_sdf.read_text().splitlines()
    counts = sdf_lines[3]
    n_atoms = int(counts[:3])
    if n_atoms != len(ligand_xyz):
        raise ValueError(
            f"Ligand atom count mismatch: SDF has {n_atoms}, PDB has {len(ligand_xyz)}. "
            "Check DOR atom ordering before running FEP."
        )
    for i, (x, y, z) in enumerate(ligand_xyz, start=4):
        old = sdf_lines[i]
        sdf_lines[i] = f"{x:10.4f}{y:10.4f}{z:10.4f}{old[30:]}"
    ligand_sdf.write_text("\n".join(sdf_lines) + "\n")
    return protein_pdb, ligand_sdf


def _phase_free_energy_kj(reporter_path: Path, temperature_k: float) -> tuple[float, float]:
    from openmm import unit
    from openmmtools.multistate import MultiStateReporter, MultiStateSamplerAnalyzer
    from openmmtools.constants import kB

    reporter = MultiStateReporter(str(reporter_path), open_mode="r")
    analyzer = MultiStateSamplerAnalyzer(reporter)
    delta_f, d_delta_f = analyzer.get_free_energy()
    kT = (kB * temperature_k * unit.kelvin).value_in_unit(unit.kilojoule_per_mole)
    dg = float(delta_f[0, -1] * kT)
    ddg = float(d_delta_f[0, -1] * kT)
    reporter.close()
    return dg, ddg


def _run_phase(htf, phase: str, output_dir: Path, args) -> Path:
    from openmm import unit
    from openmmtools import cache, mcmc, utils
    from openmmtools.multistate import MultiStateReporter
    from perses.annihilation.lambda_protocol import LambdaProtocol
    from perses.dispersed.utils import configure_platform
    from perses.samplers.multistate import HybridRepexSampler

    reporter_path = output_dir / f"{phase}.nc"
    atom_selection = args.atom_selection
    analysis_indices = htf.hybrid_topology.select(atom_selection) if atom_selection else None
    reporter = MultiStateReporter(
        str(reporter_path),
        analysis_particle_indices=analysis_indices,
        checkpoint_interval=args.checkpoint_interval,
    )
    move = mcmc.LangevinSplittingDynamicsMove(
        timestep=args.timestep_fs * unit.femtoseconds,
        collision_rate=args.collision_rate / unit.picosecond,
        n_steps=args.steps_per_cycle,
        reassign_velocities=False,
        n_restart_attempts=20,
        splitting="V R R R O R R R V",
        constraint_tolerance=1e-6,
    )
    sampler = HybridRepexSampler(
        mcmc_moves=move,
        hybrid_factory=htf,
        online_analysis_interval=args.online_analysis_interval,
    )
    sampler.setup(
        n_states=args.n_states,
        temperature=args.temperature_k * unit.kelvin,
        storage_file=reporter,
        lambda_protocol=LambdaProtocol(functions=args.lambda_protocol),
        endstates=args.sample_endstates,
    )
    platform = configure_platform(args.platform or utils.get_fastest_platform().getName())
    sampler.energy_context_cache = cache.ContextCache(capacity=None, time_to_live=None, platform=platform)
    sampler.sampler_context_cache = cache.ContextCache(capacity=None, time_to_live=None, platform=platform)
    sampler.extend(args.n_cycles)
    reporter.close()
    return reporter_path


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run exact Perses WT-to-mutant FEP cycle: ΔΔG_bind = ΔG_holo - ΔG_apo."
    )
    parser.add_argument("--mutation", required=True, help="Single mutation, e.g. V106A")
    parser.add_argument("--chain-id", default="A")
    parser.add_argument("--residue-id", default=None, help="Override residue id; defaults to parsed mutation number.")
    parser.add_argument("--old-residue", default=None, help="Override old residue three-letter code.")
    parser.add_argument("--new-residue", default=None, help="Override proposed residue three-letter code.")
    parser.add_argument("--wt-complex-pdb", type=Path, default=Path("results/md_runs/wt/rep_01/assets/wt_md_rep01_start.pdb"))
    parser.add_argument("--ligand-template-sdf", type=Path, default=Path("data/ligands/dor.sdf"))
    parser.add_argument("--ligand-resname", default="2KW")
    parser.add_argument("--output-dir", type=Path, default=Path("results/analysis/perses_point_mutation_fep"))
    parser.add_argument("--forcefield-files", default="amber14/protein.ff14SB.xml,amber14/DNA.bsc1.xml,amber14/tip3p.xml")
    parser.add_argument("--small-molecule-forcefield", default="openff-2.0.0")
    parser.add_argument("--n-states", type=int, default=11)
    parser.add_argument("--n-cycles", type=int, default=5000)
    parser.add_argument("--steps-per-cycle", type=int, default=250)
    parser.add_argument("--timestep-fs", type=float, default=4.0)
    parser.add_argument("--collision-rate", type=float, default=5.0)
    parser.add_argument("--checkpoint-interval", type=int, default=10)
    parser.add_argument("--online-analysis-interval", type=int, default=100)
    parser.add_argument("--temperature-k", type=float, default=300.0)
    parser.add_argument("--platform", default="CUDA")
    parser.add_argument("--lambda-protocol", default="default")
    parser.add_argument("--atom-selection", default="not water")
    parser.add_argument("--sample-endstates", action="store_true")
    parser.add_argument("--skip-endstate-validation", action="store_true")
    parser.add_argument("--generate-rest-htf", action="store_true")
    args = parser.parse_args()

    old_res, parsed_resid, new_res = parse_mutation_label(args.mutation)
    old_res = args.old_residue or old_res
    new_res = args.new_residue or new_res
    residue_id = args.residue_id or parsed_resid

    run_dir = args.output_dir / args.mutation / f"{args.chain_id}{residue_id}_{old_res}_to_{new_res}"
    input_dir = run_dir / "inputs"
    protein_pdb, ligand_sdf = extract_protein_and_ligand(
        args.wt_complex_pdb,
        args.ligand_template_sdf,
        input_dir,
        ligand_resname=args.ligand_resname,
    )

    start = time.time()
    from openmm import MonteCarloBarostat, unit
    from openmm import app
    from perses.app.relative_point_mutation_setup import PointMutationExecutor

    executor = PointMutationExecutor(
        protein_filename=str(protein_pdb),
        mutation_chain_id=args.chain_id,
        mutation_residue_id=str(residue_id),
        proposed_residue=new_res,
        old_residue=old_res,
        ligand_input=str(ligand_sdf),
        ligand_index=0,
        allow_undefined_stereo_sdf=True,
        is_solvated=False,
        forcefield_files=[x.strip() for x in args.forcefield_files.split(",") if x.strip()],
        small_molecule_forcefields=args.small_molecule_forcefield,
        barostat=MonteCarloBarostat(1.0 * unit.atmosphere, args.temperature_k * unit.kelvin, 50),
        forcefield_kwargs={
            "removeCMMotion": False,
            "constraints": app.HBonds,
            "hydrogenMass": 3 * unit.amus,
        },
        periodic_forcefield_kwargs={
            "nonbondedMethod": app.PME,
            "ewaldErrorTolerance": 1e-4,
        },
        conduct_endstate_validation=not args.skip_endstate_validation,
        generate_unmodified_hybrid_topology_factory=not args.generate_rest_htf,
        generate_rest_capable_hybrid_topology_factory=args.generate_rest_htf,
    )
    complex_htf = executor.get_complex_rest_htf() if args.generate_rest_htf else executor.get_complex_htf()
    apo_htf = executor.get_apo_rest_htf() if args.generate_rest_htf else executor.get_apo_htf()

    run_dir.mkdir(parents=True, exist_ok=True)
    complex_reporter = _run_phase(complex_htf, "holo", run_dir, args)
    apo_reporter = _run_phase(apo_htf, "apo", run_dir, args)
    dg_holo, ddg_holo = _phase_free_energy_kj(complex_reporter, args.temperature_k)
    dg_apo, ddg_apo = _phase_free_energy_kj(apo_reporter, args.temperature_k)
    ddg_bind = dg_holo - dg_apo
    ddg_unc = (ddg_holo**2 + ddg_apo**2) ** 0.5

    summary = {
        "mutation": args.mutation,
        "chain_id": args.chain_id,
        "residue_id": residue_id,
        "old_residue": old_res,
        "new_residue": new_res,
        "wt_complex_pdb": str(args.wt_complex_pdb),
        "protein_pdb": str(protein_pdb),
        "ligand_sdf": str(ligand_sdf),
        "dg_holo_kj_mol": dg_holo,
        "dg_holo_unc_kj_mol": ddg_holo,
        "dg_apo_kj_mol": dg_apo,
        "dg_apo_unc_kj_mol": ddg_apo,
        "ddg_bind_kj_mol": ddg_bind,
        "ddg_bind_unc_kj_mol": ddg_unc,
        "n_states": args.n_states,
        "n_cycles": args.n_cycles,
        "steps_per_cycle": args.steps_per_cycle,
        "elapsed_seconds": time.time() - start,
    }
    (run_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
