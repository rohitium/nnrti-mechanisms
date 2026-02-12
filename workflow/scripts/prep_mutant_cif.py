"""Snakemake script: apply mutation to WT CIF to produce mutant CIF."""

from pathlib import Path

from src.structure_prep.config import dor_4ncg_spec
from src.structure_prep.mutation.mutagenesis import apply_mutations
from src.structure_prep.mutation.numbering import detect_numbering_scheme
from src.structure_prep.mutation.steps import build_mutation_steps
from src.analysis.susceptibility import load_dor_susceptibilities
from src.utils import load_chain_subunits, load_residue_mappings

root = Path(".").resolve()
spec = dor_4ncg_spec(root)
mutation = snakemake.params.mutation  # noqa: F821 (injected by Snakemake)
xlsx = Path(snakemake.input.xlsx)  # noqa: F821

dor_df = load_dor_susceptibilities(xlsx, default_chain="A")
row = dor_df[dor_df["mutation"] == mutation].iloc[0]
chain_list = [c.strip().upper() for c in str(row["chain"]).split("+") if c.strip()]

chain_map = load_chain_subunits(spec.structure.cif_path)
residue_maps = load_residue_mappings(spec.structure.cif_path)
numbering = detect_numbering_scheme(spec.structure.cif_path, chain_map)

mutation_steps, _ = build_mutation_steps(
    mutation_label=mutation,
    chain_list=chain_list,
    residue_maps=residue_maps,
    numbering_scheme=numbering,
)

output_cif = Path(snakemake.output[0])  # noqa: F821
output_cif.parent.mkdir(parents=True, exist_ok=True)
apply_mutations(Path(snakemake.input.wt_cif), mutation_steps, output_cif)  # noqa: F821
