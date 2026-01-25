from __future__ import annotations

from pathlib import Path
import logging
import multiprocessing as mp
import time

import pandas as pd

from .analysis_metrics import compute_contacts, pocket_volume_proxy
from .config import dor_spec, rpv_spec
from .mutagenesis import apply_mutations
from .openmm_pipeline import compute_binding_proxy, minimize_with_restraints
from .plotting import plot_delta_metrics
from .utils import (
    ensure_dirs,
    load_chain_subunits,
    parse_mutation_group,
    project_paths,
    sanitize_label,
)


def _load_drms(drms_path: Path) -> pd.DataFrame:
    df = pd.read_csv(drms_path)
    df["drug"] = df["drug"].str.upper()
    df["mutation"] = df["mutation"].str.strip()
    if "chain" not in df.columns:
        raise ValueError("DRM table must include a 'chain' column.")
    df["chain"] = df["chain"].astype(str).str.strip().str.upper()
    df["order"] = range(len(df))
    return df


def _prepare_structure(
    cif_path: Path,
    ligand_resname: str,
    ligand_sdf: Path,
    restraint_radius: float,
    restraint_k: float,
    output_path: Path,
):
    start = time.perf_counter()
    topology, positions, forcefield = minimize_with_restraints(
        cif_path=cif_path,
        ligand_resname=ligand_resname,
        ligand_sdf=ligand_sdf,
        restraint_radius_angstrom=restraint_radius,
        restraint_k_kj_mol_nm2=restraint_k,
        output_path=output_path,
    )
    logging.info("Minimized structure in %.2fs", time.perf_counter() - start)
    start = time.perf_counter()
    energies = compute_binding_proxy(
        topology, positions, forcefield, ligand_resname=ligand_resname
    )
    logging.info("Energy proxy computed in %.2fs", time.perf_counter() - start)
    start = time.perf_counter()
    metrics_path = output_path.with_suffix(".pdb")
    contacts = compute_contacts(metrics_path, ligand_resname=ligand_resname)
    logging.info("Contacts computed in %.2fs", time.perf_counter() - start)
    start = time.perf_counter()
    pocket = pocket_volume_proxy(metrics_path, ligand_resname=ligand_resname)
    logging.info("Pocket volume computed in %.2fs", time.perf_counter() - start)
    return energies, contacts, pocket


def _mutation_worker(task: dict) -> dict:
    mutation_label = task["mutation"]
    safe_label = task["safe_label"]

    mutation_steps = parse_mutation_group(mutation_label, task["chains"])
    mut_dir = task["out_dir"] / safe_label
    ensure_dirs([mut_dir])
    mut_cif = mut_dir / f"mut_{safe_label}.cif"
    apply_mutations(
        cif_path=task["cif_path"],
        mutations=mutation_steps,
        output_path=mut_cif,
    )
    mut_min_path = mut_dir / f"mut_minimized_{safe_label}.cif"
    energies_mut, contacts_mut, pocket_mut = _prepare_structure(
        cif_path=mut_cif,
        ligand_resname=task["ligand_resname"],
        ligand_sdf=task["ligand_sdf"],
        restraint_radius=task["restraint_radius"],
        restraint_k=task["restraint_k"],
        output_path=mut_min_path,
    )

    mut_metrics = {
        "binding_proxy_kj_mol": energies_mut.binding_proxy_kj_mol,
        "contact_count": contacts_mut.contact_count,
        "hbond_count": contacts_mut.hbond_count,
        "pocket_volume_proxy": pocket_mut,
    }
    return {"base": task["base"], "mut_metrics": mut_metrics}


def _run_mutations(
    run_spec, paths, mutation_rows: pd.DataFrame, chain_map: dict[str, str]
):
    name = run_spec.structure.name
    out_dir = paths.generated / name.lower()
    ensure_dirs([out_dir, paths.results, paths.plots])

    wt_dir = out_dir / "wt"
    ensure_dirs([wt_dir])
    wt_path = wt_dir / "wt_minimized.cif"
    energies_wt, contacts_wt, pocket_wt = _prepare_structure(
        cif_path=run_spec.structure.cif_path,
        ligand_resname=run_spec.structure.ligand_resname,
        ligand_sdf=run_spec.structure.ligand_sdf,
        restraint_radius=run_spec.restraint_radius_angstrom,
        restraint_k=run_spec.restraint_k_kj_mol_nm2,
        output_path=wt_path,
    )

    wt_metrics = {
        "binding_proxy_kj_mol": energies_wt.binding_proxy_kj_mol,
        "contact_count": contacts_wt.contact_count,
        "hbond_count": contacts_wt.hbond_count,
        "pocket_volume_proxy": pocket_wt,
    }

    tasks = []
    for _, row in mutation_rows.iterrows():
        mutation_label = row["mutation"]
        chain_spec = row["chain"]
        chain_list = [c.strip().upper() for c in str(chain_spec).split("+") if c.strip()]
        if not chain_list:
            raise ValueError(f"Missing chain for mutation {mutation_label} in {name}")
        for chain_id in chain_list:
            if chain_id not in chain_map:
                raise ValueError(
                    f"Chain '{chain_id}' not found for {name} in {run_spec.structure.cif_path}"
                )
        chain_label = "+".join(chain_list)
        subunit = "+".join([chain_map.get(c, "") for c in chain_list if chain_map.get(c)])
        safe_label = sanitize_label(mutation_label)
        base = {
            "structure": name,
            "mutation": mutation_label,
            "mutation_order": int(row["order"]),
            "category": row.get("category"),
            "notes": row.get("notes"),
            "chain": chain_label,
            "subunit": subunit,
        }
        tasks.append(
            {
                "base": base,
                "mutation": mutation_label,
                "chains": chain_list,
                "safe_label": safe_label,
                "cif_path": run_spec.structure.cif_path,
                "ligand_resname": run_spec.structure.ligand_resname,
                "ligand_sdf": run_spec.structure.ligand_sdf,
                "restraint_radius": run_spec.restraint_radius_angstrom,
                "restraint_k": run_spec.restraint_k_kj_mol_nm2,
                "out_dir": out_dir,
            }
        )

    if not tasks:
        return []

    if len(tasks) == 1:
        results = [_mutation_worker(tasks[0])]
    else:
        with mp.Pool(processes=mp.cpu_count()) as pool:
            results = pool.map(_mutation_worker, tasks)

    rows = []
    for result in results:
        base = result["base"]
        mut_metrics = result["mut_metrics"]
        for state, metrics in [("WT", wt_metrics), ("MUT", mut_metrics)]:
            rows.extend(
                [
                    {
                        **base,
                        "state": state,
                        "metric": "binding_proxy_kj_mol",
                        "value": metrics["binding_proxy_kj_mol"],
                    },
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
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    root = Path(__file__).resolve().parents[1]
    paths = project_paths(root)
    ensure_dirs([paths.generated, paths.results, paths.plots])

    output_csv = paths.results / "metrics_summary.csv"
    if output_csv.exists():
        logging.info("Using existing metrics file: %s", output_csv)
        df = pd.read_csv(output_csv)
    else:
        drms = _load_drms(paths.data / "DRMs.csv")
        all_rows = []
        for spec in [rpv_spec(root), dor_spec(root)]:
            drug_rows = drms[drms["drug"] == spec.structure.name.upper()]
            if drug_rows.empty:
                logging.warning("No DRM entries found for %s", spec.structure.name)
                continue
            chain_map = load_chain_subunits(spec.structure.cif_path)
            if not chain_map:
                raise ValueError(
                    f"No RT subunit chains found in {spec.structure.cif_path}"
                )
            logging.info(
                "Chain to subunit map for %s: %s", spec.structure.name, chain_map
            )
            all_rows.extend(_run_mutations(spec, paths, drug_rows, chain_map))
        df = pd.DataFrame(all_rows)
    df.to_csv(output_csv, index=False)

    plot_delta_metrics(df, paths)


if __name__ == "__main__":
    main()
