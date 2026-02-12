"""Snakemake script: run heating + production MD for one replicate."""

import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

from src.md.manifest import MDTask
from src.md.worker import run_md_task

safe_label = snakemake.wildcards.safe_label  # noqa: F821
rep = int(snakemake.wildcards.rep)  # noqa: F821

task = MDTask(
    task_id=0,
    structure="DOR",
    mutation=snakemake.params.mutation,  # noqa: F821
    safe_label=safe_label,
    replicate=rep,
    minimized_pdb="",
    ligand_sdf=str(Path(snakemake.params.ligand_sdf).resolve()),  # noqa: F821
    ligand_resname=snakemake.params.ligand_resname,  # noqa: F821
    fold_reduction=snakemake.params.fold_reduction,  # noqa: F821
    output_json=str(Path(snakemake.output.result_json).resolve()),  # noqa: F821
    prepared_topology_pdb=str(Path(snakemake.input.topology_pdb).resolve()),  # noqa: F821
    prepared_system_xml=str(Path(snakemake.input.system_xml).resolve()),  # noqa: F821
)

run_md_task(
    task,
    heating_ps=snakemake.params.heating_ps,  # noqa: F821
    production_ns=snakemake.params.production_ns,  # noqa: F821
    report_interval=snakemake.params.report_interval,  # noqa: F821
    checkpoint_interval=snakemake.params.checkpoint_interval,  # noqa: F821
    resume_from_checkpoint=True,
)
