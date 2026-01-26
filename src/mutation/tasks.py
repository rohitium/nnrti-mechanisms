from __future__ import annotations

from itertools import combinations

import pandas as pd

from .steps import build_mutation_steps
from ..structure_prep import prepare_structure
from ..utils import ensure_dirs, sanitize_label


def compute_wt_metrics(run_spec, out_dir):
    wt_dir = out_dir / "wt"
    ensure_dirs([wt_dir])
    wt_path = wt_dir / "wt_minimized.pdb"
    energies_wt, contacts_wt, pocket_wt = prepare_structure(
        cif_path=run_spec.structure.cif_path,
        ligand_resname=run_spec.structure.ligand_resname,
        ligand_sdf=run_spec.structure.ligand_sdf,
        restraint_radius=run_spec.restraint_radius_angstrom,
        restraint_k=run_spec.restraint_k_kj_mol_nm2,
        output_path=wt_path,
    )
    return {
        "binding_proxy_kj_mol": energies_wt.binding_proxy_kj_mol,
        "contact_count": contacts_wt.contact_count,
        "hbond_count": contacts_wt.hbond_count,
        "pocket_volume_proxy": pocket_wt,
    }


def _expand_mutation_rows(mutation_rows: pd.DataFrame) -> list[dict]:
    rows = mutation_rows.to_dict("records")
    if not rows:
        return []

    def _normalize_chain_spec(chain_spec: str) -> str:
        chains = [c.strip().upper() for c in str(chain_spec).split("+") if c.strip()]
        return "+".join(chains)

    seen = {
        (row["mutation"], _normalize_chain_spec(row.get("chain", ""))) for row in rows
    }
    max_order = max(int(row["order"]) for row in rows)
    next_order = max_order + 1
    extras = []

    for row in rows:
        mutation_label = row["mutation"]
        tokens = [t.strip() for t in str(mutation_label).split("+") if t.strip()]
        if len(tokens) < 2:
            continue
        chain_spec = row.get("chain")
        chain_list = [c.strip().upper() for c in str(chain_spec).split("+") if c.strip()]
        if len(chain_list) == 1:
            chain_list = chain_list * len(tokens)
        if len(tokens) != len(chain_list):
            raise ValueError(
                f"Chain count ({len(chain_list)}) does not match mutation count ({len(tokens)}) "
                f"for {mutation_label}."
            )

        for size in range(1, len(tokens)):
            for indices in combinations(range(len(tokens)), size):
                subset_tokens = [tokens[i] for i in indices]
                subset_chains = [chain_list[i] for i in indices]
                subset_label = "+".join(subset_tokens)
                subset_chain_spec = "+".join(subset_chains)
                key = (subset_label, subset_chain_spec)
                if key in seen:
                    continue
                seen.add(key)
                derived = dict(row)
                derived["mutation"] = subset_label
                derived["chain"] = subset_chain_spec
                derived["order"] = next_order
                next_order += 1
                base_notes = derived.get("notes")
                derived_note = f"subset of {mutation_label}"
                if isinstance(base_notes, str) and base_notes.strip():
                    derived["notes"] = f"{base_notes}; {derived_note}"
                else:
                    derived["notes"] = derived_note
                extras.append(derived)

    return rows + extras


def build_tasks(
    run_spec,
    mutation_rows: pd.DataFrame,
    chain_map: dict[str, str],
    residue_maps: dict[str, dict[str, dict[str, str]]],
    numbering_scheme: dict[str, str],
    out_dir,
):
    tasks = []
    for row in _expand_mutation_rows(mutation_rows):
        mutation_label = row["mutation"]
        chain_spec = row["chain"]
        chain_list = [c.strip().upper() for c in str(chain_spec).split("+") if c.strip()]
        if not chain_list:
            raise ValueError(f"Missing chain for mutation {mutation_label}")
        for chain_id in chain_list:
            if chain_id not in chain_map:
                raise ValueError(f"Chain '{chain_id}' not found for {mutation_label}")
        subunit = "+".join(
            [chain_map.get(c, "") for c in chain_list if chain_map.get(c)]
        )
        mutation_steps, verify_targets = build_mutation_steps(
            mutation_label, chain_list, residue_maps, numbering_scheme
        )
        tasks.append(
            {
                "base": {
                    "structure": run_spec.structure.name,
                    "mutation": mutation_label,
                    "mutation_order": int(row["order"]),
                    "category": row.get("category"),
                    "notes": row.get("notes"),
                    "chain": "+".join(chain_list),
                    "subunit": subunit,
                },
                "mutation": mutation_label,
                "safe_label": sanitize_label(mutation_label),
                "mutation_steps": mutation_steps,
                "verify_targets": verify_targets,
                "cif_path": run_spec.structure.cif_path,
                "ligand_resname": run_spec.structure.ligand_resname,
                "ligand_sdf": run_spec.structure.ligand_sdf,
                "restraint_radius": run_spec.restraint_radius_angstrom,
                "restraint_k": run_spec.restraint_k_kj_mol_nm2,
                "out_dir": out_dir,
            }
        )
    return tasks
