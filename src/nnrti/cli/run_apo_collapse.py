#!/usr/bin/env python3
"""Short, densely sampled apo runs that resolve the NNIBP collapse.

The production apo trajectories save a frame every 50-500 ps, and the pocket
is already shut at the first saved frame: the apo starting structures measure
163-197 A^3, the same as the holo ones, and the first production frame is
46-80 A^3. The collapse is therefore faster than the production frame spacing
and cannot be drawn from those runs.

This restarts from the same prepared assets under the same protocol and saves
a frame every picosecond, which resolves it. Protein atoms only are written,
so a 200 ps run costs about 40 MB rather than 2 GB.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.run_apo_collapse --genotype wt --replicate 1
"""

from __future__ import annotations

import argparse
import csv
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]


def main(argv: list[str] | None = None) -> int:
    import numpy as np
    import openmm as mm
    import openmm.app as app
    import openmm.unit as u
    import MDAnalysis as mda
    from nnrti.analysis.metrics import pocket_volume_proxy_from_universe

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--genotype", required=True)
    ap.add_argument("--replicate", type=int, default=1)
    ap.add_argument("--picoseconds", type=float, default=200.0)
    ap.add_argument("--frame-ps", type=float, default=1.0)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--timestep-fs", type=float, default=2.0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--output-dir", type=Path,
                    default=ROOT / "results/analysis/apo_collapse")
    args = ap.parse_args(argv)

    g, rep = args.genotype, args.replicate
    assets = ROOT / f"results/md_runs/apo/{g}/rep_{rep:02d}/assets"
    start = assets / f"{g}_apo_md_rep{rep:02d}_start.pdb"
    sysxml = assets / f"{g}_apo_md_rep{rep:02d}_system.xml"
    for p in (start, sysxml):
        if not p.exists():
            raise SystemExit(f"missing {p}")

    pdb = app.PDBFile(str(start))
    system = mm.XmlSerializer.deserialize(sysxml.read_text())
    # same production protocol as md_protocol._create_production_simulation:
    # the barostat is added at run time, not stored in the serialised system
    system.addForce(mm.MonteCarloBarostat(1.0 * u.bar, 300 * u.kelvin, 25))
    integ = mm.LangevinMiddleIntegrator(300 * u.kelvin, 1.0 / u.picosecond,
                                        args.timestep_fs * u.femtoseconds)
    integ.setRandomNumberSeed(args.seed or (hash((g, rep)) & 0x7FFFFFFF))
    sim = app.Simulation(pdb.topology, system, integ,
                         mm.Platform.getPlatformByName("CPU"),
                         {"Threads": str(args.threads)})
    sim.context.setPositions(pdb.positions)
    sim.context.setVelocitiesToTemperature(300 * u.kelvin, args.seed or 1)

    uni = mda.Universe(str(start))
    prot = uni.select_atoms("protein").indices
    step_per_frame = max(1, int(round(args.frame_ps * 1000.0 / args.timestep_fs)))
    n_frames = int(round(args.picoseconds / args.frame_ps))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    out = args.output_dir / f"{g}_rep{rep:02d}_collapse.csv"
    t0 = time.time()
    with out.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["genotype", "replicate", "time_ps", "pocket_volume"])
        # t = 0 is the prepared structure itself, before any dynamics
        uni.atoms.positions = np.asarray(
            pdb.positions.value_in_unit(u.angstrom))
        w.writerow([g, rep, 0.0, round(pocket_volume_proxy_from_universe(uni, ""), 2)])
        fh.flush()
        for i in range(1, n_frames + 1):
            sim.step(step_per_frame)
            pos = sim.context.getState(getPositions=True).getPositions(asNumpy=True)
            uni.atoms.positions = np.asarray(pos.value_in_unit(u.angstrom))
            v = pocket_volume_proxy_from_universe(uni, "")
            w.writerow([g, rep, round(i * args.frame_ps, 3), round(float(v), 2)])
            fh.flush()
            if i % 20 == 0:
                el = time.time() - t0
                print(f"  {g} rep{rep}  {i*args.frame_ps:6.1f} ps  V={v:6.1f}  "
                      f"({el/60:.1f} min, eta {el/i*(n_frames-i)/60:.0f} min)",
                      flush=True)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
