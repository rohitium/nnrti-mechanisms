#!/usr/bin/env python3
"""Generate a Boltz affinity YAML for the WT RT + doravirine system."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.utils.cif_parser import iter_cif_loops

DEFAULT_DORAVIRINE_SMILES = (
    "CN1C(=NN=C1CN2C=CC(=NC2=O)OC3=CC(=CC=C3C#N)Cl)C(F)(F)F"
)


def _clean_sequence(raw: str) -> str:
    return re.sub(r"[^A-Za-z]", "", raw).upper()


def extract_chain_sequences(cif_path: Path) -> dict[str, str]:
    lines = cif_path.read_text().splitlines()
    chain_to_seq: dict[str, str] = {}

    for tags, data_tokens in iter_cif_loops(lines):
        if (
            "_entity_poly.pdbx_strand_id" not in tags
            or "_entity_poly.pdbx_seq_one_letter_code_can" not in tags
        ):
            continue

        strand_idx = tags.index("_entity_poly.pdbx_strand_id")
        seq_idx = tags.index("_entity_poly.pdbx_seq_one_letter_code_can")
        ncols = len(tags)

        for row in range(0, len(data_tokens), ncols):
            values = data_tokens[row : row + ncols]
            if len(values) < ncols:
                break

            seq = _clean_sequence(values[seq_idx])
            if not seq:
                continue

            for chain_id in values[strand_idx].split(","):
                chain_id = chain_id.strip()
                if not chain_id:
                    continue
                prev = chain_to_seq.get(chain_id)
                if prev is not None and prev != seq:
                    raise ValueError(
                        f"Conflicting sequences found for chain '{chain_id}' in {cif_path}"
                    )
                chain_to_seq[chain_id] = seq

    if not chain_to_seq:
        raise ValueError(
            f"Could not find polymer sequences in {cif_path} (_entity_poly loop missing)."
        )
    return dict(sorted(chain_to_seq.items(), key=lambda kv: kv[0]))


def build_yaml(
    chain_to_seq: dict[str, str],
    ligand_smiles: str,
    ligand_id: str,
    binder_id: str,
) -> str:
    if binder_id != ligand_id:
        raise ValueError("For this workflow, binder_id must match ligand_id.")

    ligand_smiles_yaml = ligand_smiles.replace("'", "''")

    lines: list[str] = ["version: 1", "sequences:"]
    for chain_id, seq in chain_to_seq.items():
        lines.extend(
            [
                "  - protein:",
                f"      id: {chain_id}",
                f"      sequence: {seq}",
            ]
        )

    lines.extend(
        [
            "  - ligand:",
            f"      id: {ligand_id}",
            f"      smiles: '{ligand_smiles_yaml}'",
            "properties:",
            "  - affinity:",
            f"      binder: {binder_id}",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-cif",
        type=Path,
        default=Path("data/prepared/dor_4ncg/wt_4ncg.cif"),
        help="WT mmCIF containing RT p66/p51 chains",
    )
    parser.add_argument(
        "--output-yaml",
        type=Path,
        default=Path("inputs/boltz/wt_rt_dor_affinity.yaml"),
        help="Output Boltz YAML path",
    )
    parser.add_argument(
        "--chains",
        default="A,B",
        help="Comma-separated chain IDs to include (default: A,B)",
    )
    parser.add_argument(
        "--ligand-id",
        default="L",
        help="Ligand ID used in the YAML sequences block",
    )
    parser.add_argument(
        "--ligand-smiles",
        default=DEFAULT_DORAVIRINE_SMILES,
        help="Doravirine SMILES string",
    )
    args = parser.parse_args()

    if not args.input_cif.exists():
        raise FileNotFoundError(f"Missing input CIF: {args.input_cif}")

    requested_chains = {tok.strip() for tok in args.chains.split(",") if tok.strip()}
    if not requested_chains:
        raise ValueError("No valid chains provided to --chains")

    chain_to_seq = extract_chain_sequences(args.input_cif)
    selected = {k: v for k, v in chain_to_seq.items() if k in requested_chains}
    missing = sorted(requested_chains.difference(selected))
    if missing:
        found = ", ".join(sorted(chain_to_seq))
        raise ValueError(
            f"Requested chains not found: {', '.join(missing)}. Available chains: {found}"
        )

    yaml_text = build_yaml(
        chain_to_seq=dict(sorted(selected.items())),
        ligand_smiles=args.ligand_smiles,
        ligand_id=args.ligand_id,
        binder_id=args.ligand_id,
    )

    args.output_yaml.parent.mkdir(parents=True, exist_ok=True)
    args.output_yaml.write_text(yaml_text)
    print(f"Wrote Boltz WT affinity input: {args.output_yaml}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
