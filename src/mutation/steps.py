from __future__ import annotations

from ..utils import one_to_three, parse_mutation_token


def build_mutation_steps(
    mutation_label: str,
    chain_list: list[str],
    residue_maps: dict[str, dict[str, dict[str, str]]],
    numbering_scheme: dict[str, str],
) -> tuple[list[tuple[str, str, str]], list[tuple[str, str, str]]]:
    tokens = mutation_label.split("+")
    if len(chain_list) == 1:
        chain_list = chain_list * len(tokens)
    if len(tokens) != len(chain_list):
        raise ValueError(
            f"Chain count ({len(chain_list)}) does not match mutation count ({len(tokens)}) "
            f"for {mutation_label}."
        )

    mutation_steps = []
    verify_targets = []
    for token, chain_id in zip(tokens, chain_list):
        old_res, resid, new_res = parse_mutation_token(token)
        expected = one_to_three(old_res)
        chain_maps = residue_maps.get(chain_id, {})
        auth_map = chain_maps.get("auth_map", {})
        if resid not in auth_map:
            raise ValueError(
                f"Residue {resid} not found in chain {chain_id} for {mutation_label}."
            )
        if auth_map[resid] != expected:
            raise ValueError(
                f"Expected {expected} at {chain_id}:{resid} for {mutation_label}, "
                f"found {auth_map[resid]}."
            )
        if numbering_scheme.get(chain_id) == "label":
            auth_to_label = chain_maps.get("auth_to_label", {})
            if resid not in auth_to_label:
                raise ValueError(
                    f"Cannot map auth residue {resid} to label_seq_id for chain {chain_id}."
                )
            resid_used = auth_to_label[resid]
        else:
            resid_used = resid
        new_three = one_to_three(new_res)
        mutation_steps.append((chain_id, resid_used, new_three))
        verify_targets.append((chain_id, expected, new_three))

    return mutation_steps, verify_targets
