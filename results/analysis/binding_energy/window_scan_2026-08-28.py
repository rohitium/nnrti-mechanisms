"""Do MM/GBSA replicate means converge as more of the trajectory is sampled?

Scores N frames in each of 5 equal time windows across the post-equilibration
region of every replicate of a few systems. If replicate means converge with
sampling, longer trajectories fix the SEM; if they stay separated, only more
replicates do.
"""
import sys, json, concurrent.futures as cf
from pathlib import Path
import numpy as np, pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
ROOT = Path("/Users/rohitpro/Career/00_Github/nnrti-mechanisms")
sys.path.insert(0, str(ROOT))
OUT = Path(__file__).resolve().parent / "win" / "window_scan.csv"

SYSTEMS = ["WT", "Y188L", "G190E", "V106A"]
N_WINDOWS = 5
PER_WINDOW = 8


def build_tasks():
    from src.analysis.result_collector import collect_md_results
    md = collect_md_results(ROOT / "manifests/md_manifest.csv", ROOT / "results")
    scr = pd.read_csv(ROOT / "results/analysis/binding_energy/frame_contact_screen.csv")
    scr = scr[(scr.status == "ok") & scr.is_clean]
    tasks = []
    for _, row in md.iterrows():
        mut, rep = str(row["mutation"]), int(row["replicate"])
        if mut not in SYSTEMS:
            continue
        clean = np.array(sorted(scr[(scr.mutation == mut) & (scr.replicate == rep)].frame.astype(int)))
        nfr = int(row.get("analysis_n_frames") or 0) or (clean.max() + 1)
        start = int(0.25 * nfr)
        clean = clean[clean >= start]
        edges = np.linspace(start, nfr, N_WINDOWS + 1).astype(int)
        for w in range(N_WINDOWS):
            band = clean[(clean >= edges[w]) & (clean < edges[w + 1])]
            if band.size == 0:
                continue
            pick = band[np.linspace(0, band.size - 1, min(PER_WINDOW, band.size)).astype(int)]
            tasks.append(dict(mutation=mut, replicate=rep, window=w,
                              frames=[int(f) for f in np.unique(pick)],
                              min_pdb=str(row["minimized_pdb"]), dcd=str(row["analysis_dcd"]),
                              topo=str(row["analysis_topology_pdb"])))
    return tasks


def run(t):
    from src.md.openmm.mmgbsa import compute_mmgbsa_from_trajectory
    try:
        mm = compute_mmgbsa_from_trajectory(
            minimized_pdb_path=Path(t["min_pdb"]), trajectory_dcd_path=Path(t["dcd"]),
            ligand_resname="2KW", ligand_sdf=ROOT / "data/ligands/dor.sdf",
            n_snapshots=len(t["frames"]), discard_fraction=0.25,
            analysis_topology_pdb_path=Path(t["topo"]), allowed_frames=t["frames"],
            snapshot_relaxation="unrestrained", relaxation_iterations=100)
        return dict(mutation=t["mutation"], replicate=t["replicate"], window=t["window"],
                    n=len(t["frames"]), binding_dg=mm.binding_dg_mean,
                    binding_dg_std=mm.binding_dg_std, vdw=mm.delta_e_vdw_mean)
    except Exception as exc:
        return dict(mutation=t["mutation"], replicate=t["replicate"], window=t["window"],
                    n=len(t["frames"]), binding_dg=np.nan, error=str(exc))


if __name__ == "__main__":
    tasks = build_tasks()
    print(f"{len(tasks)} window-tasks, {sum(len(t['frames']) for t in tasks)} frames", flush=True)
    rows = []
    with cf.ProcessPoolExecutor(max_workers=12) as ex:
        for i, r in enumerate(ex.map(run, tasks), 1):
            rows.append(r)
            print(f"[{i}/{len(tasks)}] {r['mutation']} rep{r['replicate']} w{r['window']} "
                  f"dG={r.get('binding_dg')}", flush=True)
            pd.DataFrame(rows).to_csv(OUT, index=False)
    print("wrote", OUT)
