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
    replicates: int = 1,
    jitter_seed_base: int | None = None,
    jitter_angstrom: float = 0.0,
):
    out_dir = paths.generated / run_spec.structure.name.lower()
    ensure_dirs([out_dir, paths.results, paths.plots])

    all_rows = []
    for replicate in range(1, replicates + 1):
        wt_seed = None
        if jitter_seed_base is not None:
            wt_seed = jitter_seed_base + replicate * 100000 + 1
        wt_metrics = compute_wt_metrics(
            run_spec,
            out_dir,
            replicate=replicate,
            jitter_seed=wt_seed,
            jitter_angstrom=jitter_angstrom,
        )
        tasks = build_tasks(
            run_spec,
            mutation_rows,
            chain_map,
            residue_maps,
            numbering_scheme,
            out_dir,
            replicate=replicate,
            jitter_seed_base=jitter_seed_base,
            jitter_angstrom=jitter_angstrom,
        )
        if not tasks:
            continue

        if len(tasks) == 1:
            results = [mutation_worker(tasks[0])]
        else:
            with mp.Pool(processes=mp.cpu_count()) as pool:
                results = pool.map(mutation_worker, tasks)

        all_rows.extend(build_rows(results, wt_metrics))

    return all_rows
