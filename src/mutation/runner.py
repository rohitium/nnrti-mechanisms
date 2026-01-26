from __future__ import annotations

import multiprocessing as mp

from .rows import build_rows
from .tasks import build_tasks, compute_wt_metrics
from .worker import mutation_worker
from ..utils import ensure_dirs


def run_mutations(
    run_spec,
    paths,
    mutation_rows,
    chain_map,
    residue_maps,
    numbering_scheme,
):
    out_dir = paths.generated / run_spec.structure.name.lower()
    ensure_dirs([out_dir, paths.results, paths.plots])

    wt_metrics = compute_wt_metrics(run_spec, out_dir)
    tasks = build_tasks(
        run_spec, mutation_rows, chain_map, residue_maps, numbering_scheme, out_dir
    )
    if not tasks:
        return []

    if len(tasks) == 1:
        results = [mutation_worker(tasks[0])]
    else:
        with mp.Pool(processes=mp.cpu_count()) as pool:
            results = pool.map(mutation_worker, tasks)

    return build_rows(results, wt_metrics)
