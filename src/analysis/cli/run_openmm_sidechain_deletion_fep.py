from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time

import numpy as np

from src.structure_prep.config import dor_4ncg_spec
from src.structure_prep.mutation.numbering import detect_numbering_scheme
from src.structure_prep.mutation.steps import build_mutation_steps
from src.utils import load_chain_subunits, load_residue_mappings
from src.utils.mutations import parse_mutation_token, sanitize_label


AA3_TO_1 = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
}

# Deliberately small, auditable targets.  These atoms become dummies in the WT
# topology to represent side-chain deletion endpoints.  This avoids unsafe
# arbitrary residue mapping while covering simple deletions in the panel.
SUPPORTED_DELETIONS = {
    ("ALA", "GLY"): {"CB", "HB1", "HB2", "HB3"},
    ("VAL", "ALA"): {"CG1", "HG11", "HG12", "HG13", "CG2", "HG21", "HG22", "HG23"},
    ("TYR", "PHE"): {"OH", "HH"},
}


def _protein_residues(topology, chain_id: str):
    residues = []
    for chain in topology.chains():
        if chain.id != chain_id:
            continue
        for residue in chain.residues():
            if residue.name in AA3_TO_1:
                residues.append(residue)
        return residues
    raise ValueError(f"Chain {chain_id!r} not found in topology.")


def _target_ordinal(root: Path, mutation: str, chain_id: str) -> tuple[int, str, str]:
    from openmm import app

    spec = dor_4ncg_spec(root)
    chain_map = load_chain_subunits(spec.structure.cif_path)
    residue_maps = load_residue_mappings(spec.structure.cif_path)
    numbering = detect_numbering_scheme(spec.structure.cif_path, chain_map)
    steps, _verify = build_mutation_steps(
        mutation_label=mutation,
        chain_list=[chain_id],
        residue_maps=residue_maps,
        numbering_scheme=numbering,
    )
    step_chain, step_resid, new_residue = steps[0]
    old_one, auth_resid, _new_one = parse_mutation_token(mutation)

    cif = app.PDBxFile(str(spec.structure.cif_path))
    residues = _protein_residues(cif.topology, step_chain)
    for ordinal, residue in enumerate(residues):
        if residue.id == step_resid and residue.name == residue_maps[chain_id]["auth_map"][auth_resid]:
            if AA3_TO_1[residue.name] != old_one:
                raise ValueError(
                    f"Mutation {mutation} expects {old_one} at {chain_id}:{auth_resid}, "
                    f"but cleaned CIF ordinal {ordinal} is {residue.name}."
                )
            return ordinal, residue.name, new_residue
    raise ValueError(
        f"Could not locate mutation target {mutation} in cleaned CIF chain {chain_id} "
        f"using prepared residue id {step_resid}."
    )


def _safe_label(label: str) -> str:
    return sanitize_label(label)


def _holo_mutant_pdb(label: str, replicate: int) -> Path:
    safe = _safe_label(label)
    return Path("results/md_runs") / safe / f"rep_{replicate:02d}" / "assets" / f"{safe}_md_rep{replicate:02d}_start.pdb"


def _apo_mutant_pdb(label: str, replicate: int) -> Path:
    safe = _safe_label(label).lower()
    return Path("results/md_runs/apo") / safe / f"rep_{replicate:02d}" / "assets" / f"{safe}_apo_md_rep{replicate:02d}_start.pdb"


def _wt_holo_pdb(replicate: int) -> Path:
    return Path("results/md_runs/wt") / f"rep_{replicate:02d}" / "assets" / f"wt_md_rep{replicate:02d}_start.pdb"


def _wt_holo_system(replicate: int) -> Path:
    return Path("results/md_runs/wt") / f"rep_{replicate:02d}" / "assets" / f"wt_md_rep{replicate:02d}_system.xml"


def _wt_apo_pdb(replicate: int) -> Path:
    return Path("results/md_runs/apo/wt") / f"rep_{replicate:02d}" / "assets" / f"wt_apo_md_rep{replicate:02d}_start.pdb"


def _wt_apo_system(replicate: int) -> Path:
    return Path("results/md_runs/apo/wt") / f"rep_{replicate:02d}" / "assets" / f"wt_apo_md_rep{replicate:02d}_system.xml"


def _classify_mutation(root: Path, mutation: str, chain_id: str, replicate: int) -> dict:
    tokens = [token.strip().upper() for token in mutation.split("+") if token.strip()]
    row = {
        "mutation": mutation,
        "n_sites": len(tokens),
        "supported": False,
        "reason": "",
        "old_residue": "",
        "new_residue": "",
        "cleaned_cif_residue_ordinal": "",
    }
    if len(tokens) != 1:
        row["reason"] = "multi-site mutations are not yet supported by the custom deletion-only driver"
        return row
    try:
        ordinal, old_residue, new_residue = _target_ordinal(root, mutation, chain_id)
    except Exception as exc:
        row["reason"] = f"mapping failed: {exc}"
        return row
    row.update(
        {
            "old_residue": old_residue,
            "new_residue": new_residue,
            "cleaned_cif_residue_ordinal": ordinal,
        }
    )
    if (old_residue, new_residue) not in SUPPORTED_DELETIONS:
        row["reason"] = f"{old_residue}->{new_residue} is not a supported side-chain deletion"
        return row
    for path in (_holo_mutant_pdb(mutation, replicate), _apo_mutant_pdb(mutation, replicate)):
        if not path.exists():
            row["reason"] = f"missing prepared endpoint asset: {path}"
            return row
        try:
            _validate_prepared_residue(path, chain_id, ordinal, new_residue)
        except Exception as exc:
            row["reason"] = f"prepared endpoint validation failed for {path}: {exc}"
            return row
    row["supported"] = True
    row["reason"] = "supported side-chain deletion with validated prepared endpoints"
    return row


def _read_topology(path: Path):
    from openmm import app

    suffix = path.suffix.lower()
    if suffix == ".cif":
        return app.PDBxFile(str(path))
    return app.PDBFile(str(path))


def _validate_prepared_residue(
    pdb_path: Path,
    chain_id: str,
    ordinal: int,
    expected_residue: str,
) -> tuple[object, object]:
    pdb = _read_topology(pdb_path)
    residues = _protein_residues(pdb.topology, chain_id)
    if ordinal >= len(residues):
        raise ValueError(f"{pdb_path} chain {chain_id} has only {len(residues)} protein residues.")
    residue = residues[ordinal]
    if residue.name != expected_residue:
        raise ValueError(
            f"{pdb_path} chain {chain_id} ordinal {ordinal} is {residue.name} residue id {residue.id}; "
            f"expected {expected_residue}."
        )
    return pdb, residue


def _alchemical_atom_indices(residue, old_residue: str, new_residue: str) -> list[int]:
    delete_names = SUPPORTED_DELETIONS.get((old_residue, new_residue))
    if delete_names is None:
        supported = ", ".join(f"{a}->{b}" for a, b in SUPPORTED_DELETIONS)
        raise ValueError(
            f"Custom OpenMM side-chain deletion FEP currently supports only: {supported}. "
            f"Requested {old_residue}->{new_residue}."
        )
    atoms = {atom.name: atom.index for atom in residue.atoms()}
    missing = sorted(delete_names.difference(atoms))
    if missing:
        raise ValueError(f"Prepared WT residue {residue} lacks expected alchemical atoms: {missing}")
    return [atoms[name] for name in sorted(delete_names)]


def _deserialize_system(path: Path):
    from openmm import XmlSerializer

    return XmlSerializer.deserialize(path.read_text())


def _find_nonbonded_force(system):
    from openmm import NonbondedForce

    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            return force
    raise ValueError("System does not contain an OpenMM NonbondedForce.")


def _capture_nonbonded_parameters(nonbonded):
    particles = [nonbonded.getParticleParameters(i) for i in range(nonbonded.getNumParticles())]
    exceptions = [nonbonded.getExceptionParameters(i) for i in range(nonbonded.getNumExceptions())]
    return particles, exceptions


def _set_alchemical_nonbonded(nonbonded, original, alchemical_atoms: set[int], strength: float) -> None:
    particles, exceptions = original
    charge_scale = math.sqrt(strength)
    lj_scale = strength
    for atom_index, (charge, sigma, epsilon) in enumerate(particles):
        if atom_index in alchemical_atoms:
            nonbonded.setParticleParameters(atom_index, charge * charge_scale, sigma, epsilon * lj_scale)
        else:
            nonbonded.setParticleParameters(atom_index, charge, sigma, epsilon)
    for exception_index, (i, j, charge_prod, sigma, epsilon) in enumerate(exceptions):
        if int(i) in alchemical_atoms or int(j) in alchemical_atoms:
            nonbonded.setExceptionParameters(
                exception_index,
                i,
                j,
                charge_prod * strength,
                sigma,
                epsilon * lj_scale,
            )
        else:
            nonbonded.setExceptionParameters(exception_index, i, j, charge_prod, sigma, epsilon)


def _potential_kj(context) -> float:
    from openmm import unit

    state = context.getState(getEnergy=True)
    return float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))


def _make_context(system, positions, platform_name: str, temperature_k: float, seed: int):
    from openmm import LangevinMiddleIntegrator, Platform, unit
    from openmm import app

    integrator = LangevinMiddleIntegrator(
        temperature_k * unit.kelvin,
        1.0 / unit.picosecond,
        2.0 * unit.femtoseconds,
    )
    integrator.setRandomNumberSeed(seed)
    platform = Platform.getPlatformByName(platform_name)
    simulation = app.Simulation(None, system, integrator, platform)
    simulation.context.setPositions(positions)
    return simulation


def _run_phase(
    phase: str,
    pdb_path: Path,
    system_path: Path,
    chain_id: str,
    ordinal: int,
    old_residue: str,
    new_residue: str,
    args,
) -> dict:
    from openmm import unit

    pdb, residue = _validate_prepared_residue(pdb_path, chain_id, ordinal, old_residue)
    alchemical_indices = _alchemical_atom_indices(residue, old_residue, new_residue)
    lambda_values = np.linspace(1.0, 0.0, args.n_windows + 1)
    beta = 1.0 / (unit.MOLAR_GAS_CONSTANT_R * args.temperature_k * unit.kelvin).value_in_unit(
        unit.kilojoule_per_mole
    )

    windows = []
    total_dg = 0.0
    for window_index, (lam_a, lam_b) in enumerate(zip(lambda_values[:-1], lambda_values[1:])):
        system = _deserialize_system(system_path)
        nonbonded = _find_nonbonded_force(system)
        original = _capture_nonbonded_parameters(nonbonded)
        atoms = set(alchemical_indices)
        _set_alchemical_nonbonded(nonbonded, original, atoms, float(lam_a))
        simulation = _make_context(
            system=system,
            positions=pdb.positions,
            platform_name=args.platform,
            temperature_k=args.temperature_k,
            seed=args.seed + window_index,
        )
        context = simulation.context
        if args.equilibration_steps:
            simulation.step(args.equilibration_steps)

        delta_u = []
        for _sample in range(args.samples_per_window):
            if args.steps_per_sample:
                simulation.step(args.steps_per_sample)
            u_a = _potential_kj(context)
            _set_alchemical_nonbonded(nonbonded, original, atoms, float(lam_b))
            nonbonded.updateParametersInContext(context)
            u_b = _potential_kj(context)
            _set_alchemical_nonbonded(nonbonded, original, atoms, float(lam_a))
            nonbonded.updateParametersInContext(context)
            delta_u.append(beta * (u_b - u_a))

        delta_u_arr = np.asarray(delta_u, dtype=float)
        max_arg = float(np.max(-delta_u_arr))
        mean_exp = float(np.exp(-delta_u_arr - max_arg).mean() * math.exp(max_arg))
        dg = -math.log(mean_exp) / beta
        sem_du = float(delta_u_arr.std(ddof=1) / math.sqrt(len(delta_u_arr))) if len(delta_u_arr) > 1 else 0.0
        total_dg += dg
        windows.append(
            {
                "phase": phase,
                "window": window_index,
                "lambda_a": float(lam_a),
                "lambda_b": float(lam_b),
                "n_samples": int(len(delta_u_arr)),
                "dg_kj_mol": float(dg),
                "delta_u_mean": float(delta_u_arr.mean()),
                "delta_u_sem": sem_du,
            }
        )
    return {
        "phase": phase,
        "pdb_path": str(pdb_path),
        "system_path": str(system_path),
        "prepared_residue_id": residue.id,
        "prepared_residue_name": residue.name,
        "alchemical_atom_indices": alchemical_indices,
        "alchemical_atom_names": [atom.name for atom in residue.atoms() if atom.index in set(alchemical_indices)],
        "dg_kj_mol": float(total_dg),
        "windows": windows,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Run a minimal OpenMM side-chain deletion FEP thermodynamic cycle. "
            "Only auditable deletion endpoints with validated prepared WT/mutant assets are run."
        )
    )
    parser.add_argument("--mutation", default="V106A", help="Mutation label, comma-list, or 'all'.")
    parser.add_argument("--chain-id", default="A")
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--susceptibility-xlsx", type=Path, default=Path("data/DRM-susceptibilities.csv.xlsx"))
    parser.add_argument("--replicate", type=int, default=1)
    parser.add_argument("--holo-wt-pdb", type=Path, default=None)
    parser.add_argument("--holo-wt-system", type=Path, default=None)
    parser.add_argument("--apo-wt-pdb", type=Path, default=None)
    parser.add_argument("--apo-wt-system", type=Path, default=None)
    parser.add_argument("--classify-only", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=Path("results/analysis/openmm_sidechain_deletion_fep"))
    parser.add_argument("--n-windows", type=int, default=4)
    parser.add_argument("--equilibration-steps", type=int, default=100)
    parser.add_argument("--samples-per-window", type=int, default=10)
    parser.add_argument("--steps-per-sample", type=int, default=50)
    parser.add_argument("--temperature-k", type=float, default=300.0)
    parser.add_argument("--platform", default="CPU")
    parser.add_argument("--seed", type=int, default=20260720)
    args = parser.parse_args()

    root = args.root.resolve()
    if args.mutation.strip().lower() == "all":
        from src.analysis.susceptibility import load_dor_susceptibilities

        panel = load_dor_susceptibilities(args.susceptibility_xlsx, default_chain=args.chain_id)
        mutations = panel["mutation"].astype(str).tolist()
    else:
        mutations = [token.strip() for token in args.mutation.split(",") if token.strip()]

    args.holo_wt_pdb = args.holo_wt_pdb or _wt_holo_pdb(args.replicate)
    args.holo_wt_system = args.holo_wt_system or _wt_holo_system(args.replicate)
    args.apo_wt_pdb = args.apo_wt_pdb or _wt_apo_pdb(args.replicate)
    args.apo_wt_system = args.apo_wt_system or _wt_apo_system(args.replicate)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    support_rows = [_classify_mutation(root, mutation, args.chain_id, args.replicate) for mutation in mutations]
    try:
        import pandas as pd

        pd.DataFrame(support_rows).to_csv(args.output_dir / "mutation_support.csv", index=False)
    except Exception:
        (args.output_dir / "mutation_support.json").write_text(json.dumps(support_rows, indent=2) + "\n")

    if args.classify_only:
        print(json.dumps({"n_mutations": len(support_rows), "supported": sum(bool(r["supported"]) for r in support_rows)}, indent=2))
        print(f"Wrote {args.output_dir / 'mutation_support.csv'}")
        return 0

    summaries = []
    for support in support_rows:
        if not support["supported"]:
            print(f"SKIP {support['mutation']}: {support['reason']}")
            continue
        mutation = str(support["mutation"])
        ordinal = int(support["cleaned_cif_residue_ordinal"])
        old_residue = str(support["old_residue"])
        new_residue = str(support["new_residue"])

        start = time.time()
        holo = _run_phase(
            "holo",
            args.holo_wt_pdb,
            args.holo_wt_system,
            args.chain_id,
            ordinal,
            old_residue,
            new_residue,
            args,
        )
        apo = _run_phase(
            "apo",
            args.apo_wt_pdb,
            args.apo_wt_system,
            args.chain_id,
            ordinal,
            old_residue,
            new_residue,
            args,
        )
        ddg = holo["dg_kj_mol"] - apo["dg_kj_mol"]

        out_dir = args.output_dir / mutation
        out_dir.mkdir(parents=True, exist_ok=True)
        summary = {
            "mutation": mutation,
            "chain_id": args.chain_id,
            "cleaned_cif_residue_ordinal": ordinal,
            "old_residue": old_residue,
            "new_residue": new_residue,
            "method": "single-topology side-chain deletion FEP; NonbondedForce charges/LJ for deleted atoms scaled to dummy endpoint",
            "dg_holo_kj_mol": holo["dg_kj_mol"],
            "dg_apo_kj_mol": apo["dg_kj_mol"],
            "ddg_bind_kj_mol": float(ddg),
            "ddg_bind_kcal_mol": float(ddg / 4.184),
            "elapsed_seconds": time.time() - start,
            "settings": {
                "replicate": args.replicate,
                "n_windows": args.n_windows,
                "equilibration_steps": args.equilibration_steps,
                "samples_per_window": args.samples_per_window,
                "steps_per_sample": args.steps_per_sample,
                "temperature_k": args.temperature_k,
                "platform": args.platform,
                "seed": args.seed,
            },
            "holo": holo,
            "apo": apo,
        }
        (out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        rows = holo["windows"] + apo["windows"]
        try:
            import pandas as pd

            pd.DataFrame(rows).to_csv(out_dir / "windows.csv", index=False)
        except Exception:
            pass
        summaries.append(summary)
        print(json.dumps({k: summary[k] for k in ["mutation", "dg_holo_kj_mol", "dg_apo_kj_mol", "ddg_bind_kj_mol", "ddg_bind_kcal_mol"]}, indent=2))
        print(f"Wrote {out_dir / 'summary.json'}")

    if summaries:
        try:
            import pandas as pd

            pd.DataFrame(
                [
                    {
                        "mutation": s["mutation"],
                        "dg_holo_kj_mol": s["dg_holo_kj_mol"],
                        "dg_apo_kj_mol": s["dg_apo_kj_mol"],
                        "ddg_bind_kj_mol": s["ddg_bind_kj_mol"],
                        "ddg_bind_kcal_mol": s["ddg_bind_kcal_mol"],
                    }
                    for s in summaries
                ]
            ).to_csv(args.output_dir / "summary.csv", index=False)
        except Exception:
            pass
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
