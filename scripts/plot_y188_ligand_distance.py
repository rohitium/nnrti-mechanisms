#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _choose_target_residue_near_ligand(
    universe,
    ligand_sel: str,
    target_resid: int,
    target_segid: str | None,
    expected_resname: str | None,
):
    from MDAnalysis.lib.distances import distance_array

    ligand = universe.select_atoms(ligand_sel)
    if ligand.n_atoms == 0:
        raise ValueError(f"Ligand selection returned 0 atoms: {ligand_sel}")

    sel_parts = ["protein", f"resid {int(target_resid)}"]
    if target_segid:
        sel_parts.append(f"segid {target_segid}")
    if expected_resname:
        sel_parts.append(f"resname {expected_resname}")
    candidates = universe.select_atoms(" and ".join(sel_parts) + " and not name H*")
    residues = list(candidates.residues)
    if not residues:
        raise ValueError(
            f"No protein residue found for resid={target_resid}, segid={target_segid}, resname={expected_resname}"
        )

    universe.trajectory[0]
    best = None
    best_d = np.inf
    for residue in residues:
        heavy = residue.atoms.select_atoms("not name H*")
        if heavy.n_atoms == 0:
            continue
        d = distance_array(heavy.positions, ligand.positions).min()
        if d < best_d:
            best_d = float(d)
            best = residue
    if best is None:
        raise ValueError("Could not find heavy atoms for resid 188 candidates.")
    return best


def _distance_timeseries(
    topology: Path,
    trajectory: Path,
    label: str,
    ligand_resname: str,
    frame_stride: int,
    target_resid: int,
    target_segid: str | None,
    expected_resname: str | None,
):
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import distance_array

    u = mda.Universe(str(topology), str(trajectory))
    lig = u.select_atoms(f"resname {ligand_resname} and not name H*")
    if lig.n_atoms == 0:
        raise ValueError(f"No ligand atoms found with resname {ligand_resname}")

    target_res = _choose_target_residue_near_ligand(
        u,
        f"resname {ligand_resname} and not name H*",
        target_resid=target_resid,
        target_segid=target_segid,
        expected_resname=expected_resname,
    )
    target_sidechain = target_res.atoms.select_atoms("not name N CA C O OXT and not name H*")
    if target_sidechain.n_atoms == 0:
        target_sidechain = target_res.atoms.select_atoms("not name H*")
    if target_sidechain.n_atoms == 0:
        raise ValueError("Target residue has no heavy atoms.")

    rows: list[dict] = []
    for ts in u.trajectory[:: max(1, frame_stride)]:
        dmin = float(distance_array(target_sidechain.positions, lig.positions).min())
        time_ps = float(ts.time) if ts.time is not None else float(ts.frame)
        rows.append(
            {
                "system": label,
                "frame": int(ts.frame),
                "time_ps": time_ps,
                "time_ns": time_ps / 1000.0,
                "target_resname": str(target_res.resname),
                "target_resid": int(target_res.resid),
                "target_segid": str(target_res.segid),
                "n_target_atoms": int(target_sidechain.n_atoms),
                "ligand_resname": ligand_resname,
                "min_sidechain_ligand_distance_angstrom": dmin,
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compute and plot Y188/L188 sidechain-to-ligand distance trajectories."
    )
    parser.add_argument("--wt-topology", type=Path, required=True)
    parser.add_argument("--wt-trajectory", type=Path, required=True)
    parser.add_argument("--mut-topology", type=Path, required=True)
    parser.add_argument("--mut-trajectory", type=Path, required=True)
    parser.add_argument("--ligand-resname", type=str, default="2KW")
    parser.add_argument("--target-resid", type=int, default=188)
    parser.add_argument("--target-segid", type=str, default=None)
    parser.add_argument("--wt-resname", type=str, default=None)
    parser.add_argument("--mut-resname", type=str, default=None)
    parser.add_argument("--mut-label", type=str, default="Y188L")
    parser.add_argument("--frame-stride", type=int, default=1)
    parser.add_argument("--out-csv", type=Path, required=True)
    parser.add_argument("--out-plot", type=Path, required=True)
    args = parser.parse_args()

    wt_df = _distance_timeseries(
        topology=args.wt_topology,
        trajectory=args.wt_trajectory,
        label="WT",
        ligand_resname=args.ligand_resname,
        frame_stride=args.frame_stride,
        target_resid=args.target_resid,
        target_segid=args.target_segid,
        expected_resname=args.wt_resname,
    )
    mut_df = _distance_timeseries(
        topology=args.mut_topology,
        trajectory=args.mut_trajectory,
        label=args.mut_label,
        ligand_resname=args.ligand_resname,
        frame_stride=args.frame_stride,
        target_resid=args.target_resid,
        target_segid=args.target_segid,
        expected_resname=args.mut_resname,
    )
    df = pd.concat([wt_df, mut_df], ignore_index=True)

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    args.out_plot.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out_csv, index=False)

    plt.figure(figsize=(10, 4.8))
    for name, group in df.groupby("system"):
        g = group.sort_values("time_ns")
        mean_val = float(g["min_sidechain_ligand_distance_angstrom"].mean())
        (line,) = plt.plot(
            g["time_ns"].to_numpy(),
            g["min_sidechain_ligand_distance_angstrom"].to_numpy(),
            linewidth=1.5,
            label=name,
        )
        plt.axhline(
            y=mean_val,
            linestyle="--",
            linewidth=1.2,
            alpha=0.8,
            color=line.get_color(),
            label=f"{name} mean ({mean_val:.2f} A)",
        )
    plt.xlabel("Time (ns)")
    plt.ylabel("Min distance: Y188/Y188L sidechain to 2KW (A)")
    plt.title("Y188/Y188L site sidechain-ligand distance trajectory")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(args.out_plot, dpi=300)

    print(f"Wrote {args.out_csv}")
    print(f"Wrote {args.out_plot}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
