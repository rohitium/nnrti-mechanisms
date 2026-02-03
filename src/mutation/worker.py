from __future__ import annotations

from .mutagenesis import apply_mutations
from ..structure_prep import compute_alchemical_binding_metric, prepare_structure
from .verify import verify_mutations
from ..utils import ensure_dirs


def mutation_worker(task: dict) -> dict:
    mutation_label = task["mutation"]
    safe_label = task["safe_label"]
    mutation_steps = task["mutation_steps"]
    replicate = task.get("replicate", 1)

    mut_dir = task["out_dir"] / safe_label / f"rep_{replicate:02d}"
    ensure_dirs([mut_dir])
    mut_cif = mut_dir / f"mut_{safe_label}_rep{replicate:02d}.cif"
    apply_mutations(
        cif_path=task["cif_path"],
        mutations=mutation_steps,
        output_path=mut_cif,
    )
    verify_mutations(task["cif_path"], mut_cif, task["verify_targets"])
    mut_min_path = mut_dir / f"mut_minimized_{safe_label}_rep{replicate:02d}.pdb"
    energies_mut, contacts_mut, pocket_mut = prepare_structure(
        cif_path=mut_cif,
        ligand_resname=task["ligand_resname"],
        ligand_sdf=task["ligand_sdf"],
        restraint_radius=task["restraint_radius"],
        restraint_k=task["restraint_k"],
        output_path=mut_min_path,
        jitter_seed=task.get("jitter_seed"),
        jitter_angstrom=task.get("jitter_angstrom", 0.0),
    )
    mut_alchemy_json = mut_dir / f"mut_alchemy_{safe_label}_rep{replicate:02d}.json"
    mut_binding_dg = compute_alchemical_binding_metric(
        minimized_pdb_path=mut_min_path,
        ligand_resname=task["ligand_resname"],
        ligand_sdf=task["ligand_sdf"],
        config=task["alchemy_config"],
        output_json=mut_alchemy_json,
        metadata={
            **task["base"],
            "state": "MUT",
        },
    )

    mut_metrics = {
        "binding_delta_g_kj_mol": mut_binding_dg,
        "binding_proxy_kj_mol": energies_mut.binding_proxy_kj_mol,
        "contact_count": contacts_mut.contact_count,
        "hbond_count": contacts_mut.hbond_count,
        "pocket_volume_proxy": pocket_mut,
    }
    return {"base": task["base"], "mut_metrics": mut_metrics}
