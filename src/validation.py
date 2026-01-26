from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from .config import dor_spec, rpv_spec
from .numbering import detect_numbering_scheme
from .mutation.mutagenesis import apply_mutations
from .mutation.steps import build_mutation_steps
from .mutation.verify import verify_mutations
from .utils import load_chain_subunits, load_residue_mappings, one_to_three, parse_mutation_token


def validate_mutations(drms: pd.DataFrame, root: Path) -> None:
    for spec in [rpv_spec(root), dor_spec(root)]:
        drug_rows = drms[drms["drug"] == spec.structure.name.upper()]
        if drug_rows.empty:
            logging.warning("No DRM entries found for %s", spec.structure.name)
            continue
        chain_map = load_chain_subunits(spec.structure.cif_path)
        residue_maps = load_residue_mappings(spec.structure.cif_path)
        numbering_scheme = detect_numbering_scheme(spec.structure.cif_path, chain_map)

        total = 0
        for _, row in drug_rows.iterrows():
            mutation_label = row["mutation"]
            chain_spec = row["chain"]
            chain_list = [
                c.strip().upper() for c in str(chain_spec).split("+") if c.strip()
            ]
            tokens = mutation_label.split("+")
            if len(chain_list) == 1:
                chain_list = chain_list * len(tokens)
            if len(tokens) != len(chain_list):
                raise ValueError(
                    f"Chain count ({len(chain_list)}) does not match mutation count ({len(tokens)}) "
                    f"for {mutation_label}."
                )
            for token, chain_id in zip(tokens, chain_list):
                old_res, resid, new_res = parse_mutation_token(token)
                expected = one_to_three(old_res)
                chain_maps = residue_maps.get(chain_id, {})
                auth_map = chain_maps.get("auth_map", {})
                if resid not in auth_map:
                    raise ValueError(
                        f"Residue {resid} not found in chain {chain_id} for {mutation_label}."
                    )
                found = auth_map[resid]
                if found != expected:
                    raise ValueError(
                        f"Expected {expected} at {chain_id}:{resid} for {mutation_label}, "
                        f"found {found}."
                    )
                resid_used = resid
                if numbering_scheme.get(chain_id) == "label":
                    auth_to_label = chain_maps.get("auth_to_label", {})
                    if resid not in auth_to_label:
                        raise ValueError(
                            f"Cannot map auth residue {resid} to label_seq_id for chain {chain_id}."
                        )
                    resid_used = auth_to_label[resid]
                total += 1
                logging.info(
                    "Validated %s %s: chain %s auth %s %s->%s (apply id %s; scheme %s)",
                    spec.structure.name,
                    mutation_label,
                    chain_id,
                    resid,
                    expected,
                    one_to_three(new_res),
                    resid_used,
                    numbering_scheme.get(chain_id),
                )
        logging.info("Validated %d substitutions for %s", total, spec.structure.name)
    logging.info(
        "Validation OK: all DRM substitutions are consistent with the CIF sequences."
    )


def verify_mutations_only(drms: pd.DataFrame, root: Path) -> None:
    for spec in [rpv_spec(root), dor_spec(root)]:
        drug_rows = drms[drms["drug"] == spec.structure.name.upper()]
        if drug_rows.empty:
            logging.warning("No DRM entries found for %s", spec.structure.name)
            continue
        chain_map = load_chain_subunits(spec.structure.cif_path)
        residue_maps = load_residue_mappings(spec.structure.cif_path)
        numbering_scheme = detect_numbering_scheme(spec.structure.cif_path, chain_map)

        for _, row in drug_rows.iterrows():
            mutation_label = row["mutation"]
            chain_spec = row["chain"]
            chain_list = [
                c.strip().upper() for c in str(chain_spec).split("+") if c.strip()
            ]
            steps, targets = build_mutation_steps(
                mutation_label, chain_list, residue_maps, numbering_scheme
            )
            tmp_path = (
                root
                / "data"
                / "generated"
                / "_verify"
                / spec.structure.name.lower()
                / f"{mutation_label}.cif"
            )
            tmp_path.parent.mkdir(parents=True, exist_ok=True)
            apply_mutations(spec.structure.cif_path, steps, tmp_path)
            verify_mutations(spec.structure.cif_path, tmp_path, targets)
            logging.info("Verified %s %s", spec.structure.name, mutation_label)
    logging.info("Verification OK: all DRM substitutions are applied correctly.")
