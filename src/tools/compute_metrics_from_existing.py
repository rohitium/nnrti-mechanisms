from __future__ import annotations

import argparse
import logging
from pathlib import Path
import sys
import time

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.analysis_metrics import compute_contacts, pocket_volume_proxy
from src.config import dor_spec, rpv_spec
from src.drm_io import load_drms
from src.mutation.tasks import build_tasks
from src.numbering import detect_numbering_scheme
from src.openmm.energy import compute_binding_proxy
from src.openmm.ligand import build_forcefield, load_ligand_molecule
from src.openmm.require import require_module
from src.utils import load_chain_subunits, load_residue_mappings, project_paths


def _compute_metrics_from_pdb(
    pdb_path: Path,
    ligand_resname: str,
    ligand_sdf: Path,
    skip_energy: bool,
) -> dict:
    logging.info("Metrics: %s", pdb_path)
    metrics: dict[str, float] = {}
    if not skip_energy:
        start = time.perf_counter()
        app = require_module("openmm.app")
        with open(pdb_path, "r") as handle:
            pdb = app.PDBFile(handle)
        ligand = load_ligand_molecule(ligand_sdf)
        forcefield = build_forcefield([ligand])
        energies = compute_binding_proxy(
            pdb.topology, pdb.positions, forcefield, ligand_resname=ligand_resname
        )
        metrics["binding_proxy_kj_mol"] = energies.binding_proxy_kj_mol
        logging.info("Energy computed in %.2fs", time.perf_counter() - start)

    start = time.perf_counter()
    contacts = compute_contacts(pdb_path, ligand_resname=ligand_resname)
    metrics["contact_count"] = contacts.contact_count
    metrics["hbond_count"] = contacts.hbond_count
    logging.info("Contacts/H-bonds in %.2fs", time.perf_counter() - start)
    start = time.perf_counter()
    metrics["pocket_volume_proxy"] = pocket_volume_proxy(
        pdb_path, ligand_resname=ligand_resname
    )
    logging.info("Pocket volume in %.2fs", time.perf_counter() - start)
    return metrics


def _rows_from_existing(
    run_spec,
    paths,
    mutation_rows,
    chain_map,
    residue_maps,
    numbering_scheme,
    replicates: int,
    skip_energy: bool,
):
    out_dir = paths.generated / run_spec.structure.name.lower()
    rows: list[dict] = []
    logging.info("Structure: %s", run_spec.structure.name)
    for rep in range(1, replicates + 1):
        logging.info("Replicate %d/%d", rep, replicates)
        wt_path = (
            out_dir / "wt" / f"rep_{rep:02d}" / f"wt_minimized_rep{rep:02d}.pdb"
        )
        if not wt_path.exists():
            logging.warning("Missing WT: %s", wt_path)
            continue
        wt_metrics = _compute_metrics_from_pdb(
            wt_path, run_spec.structure.ligand_resname, run_spec.structure.ligand_sdf, skip_energy
        )
        tasks = build_tasks(
            run_spec,
            mutation_rows,
            chain_map,
            residue_maps,
            numbering_scheme,
            out_dir,
            replicate=rep,
            jitter_seed_base=None,
            jitter_angstrom=0.0,
        )
        logging.info("Tasks found: %d", len(tasks))
        for task in tasks:
            safe_label = task["safe_label"]
            mut_path = (
                out_dir
                / safe_label
                / f"rep_{rep:02d}"
                / f"mut_minimized_{safe_label}_rep{rep:02d}.pdb"
            )
            if not mut_path.exists():
                logging.info("Skip (missing): %s", mut_path)
                continue
            logging.info("Mutation: %s", task["mutation"])
            mut_metrics = _compute_metrics_from_pdb(
                mut_path,
                run_spec.structure.ligand_resname,
                run_spec.structure.ligand_sdf,
                skip_energy,
            )
            base = task["base"]
            for state, metrics in [("WT", wt_metrics), ("MUT", mut_metrics)]:
                if "binding_proxy_kj_mol" in metrics:
                    rows.append(
                        {
                            **base,
                            "state": state,
                            "metric": "binding_proxy_kj_mol",
                            "value": metrics["binding_proxy_kj_mol"],
                        }
                    )
                rows.extend(
                    [
                        {
                            **base,
                            "state": state,
                            "metric": "contact_count",
                            "value": metrics["contact_count"],
                        },
                        {
                            **base,
                            "state": state,
                            "metric": "hbond_count",
                            "value": metrics["hbond_count"]
                            if metrics["hbond_count"] is not None
                            else float("nan"),
                        },
                        {
                            **base,
                            "state": state,
                            "metric": "pocket_volume_proxy",
                            "value": metrics["pocket_volume_proxy"],
                        },
                    ]
                )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute metrics from already minimized structures."
    )
    parser.add_argument("--replicates", type=int, default=3)
    parser.add_argument(
        "--skip-energy",
        action="store_true",
        help="Skip binding proxy energy calculation (faster).",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output CSV path (default: results/metrics_summary.csv).",
    )
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )

    root = Path.cwd()
    paths = project_paths(root)

    drms = load_drms(paths.data / "DRMs.csv")
    output_csv = Path(args.output) if args.output else paths.results / "metrics_summary.csv"
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    wrote_header = output_csv.exists()
    for spec in [rpv_spec(root), dor_spec(root)]:
        drug_rows = drms[drms["drug"] == spec.structure.name.upper()]
        if drug_rows.empty:
            continue
        chain_map = load_chain_subunits(spec.structure.cif_path)
        residue_maps = load_residue_mappings(spec.structure.cif_path)
        numbering = detect_numbering_scheme(spec.structure.cif_path, chain_map)
        rows = _rows_from_existing(
            spec,
            paths,
            drug_rows,
            chain_map,
            residue_maps,
            numbering,
            replicates=args.replicates,
            skip_energy=args.skip_energy,
        )
        if not rows:
            continue
        df = pd.DataFrame(rows)
        df.to_csv(output_csv, mode="a", header=not wrote_header, index=False)
        wrote_header = True
        logging.info("Appended %d rows to %s", len(df), output_csv)

    print(f"Done. Output at {output_csv}")


if __name__ == "__main__":
    main()
