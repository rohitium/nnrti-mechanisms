#!/usr/bin/env python3
"""Short, densely sampled apo run: save protein coordinates every picosecond.

Deliberately depends on OpenMM and numpy only, so it runs under Sherlock's
py-openmm module without building an environment. No pocket volume is computed
here: the coordinates come back and are measured with the same routine used
for every other number in the paper, so the comparison cannot drift.

Protocol matches md_protocol._create_production_simulation: the barostat is
added at run time rather than stored in the serialised system.

Usage
-----
    python collapse_run.py GENOTYPE REPLICATE ASSET_ROOT OUT_DIR [PICOSECONDS]
"""
import sys, time
from pathlib import Path

import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as u

SOLVENT = {"HOH", "WAT", "TIP3", "NA", "CL", "K", "MG", "SOD", "CLA"}


def main() -> int:
    g, rep = sys.argv[1], int(sys.argv[2])
    assets = Path(sys.argv[3]) / g / f"rep_{rep:02d}" / "assets"
    out = Path(sys.argv[4])
    ps_total = float(sys.argv[5]) if len(sys.argv) > 5 else 200.0
    frame_ps, dt_fs = 1.0, 2.0

    start = assets / f"{g}_apo_md_rep{rep:02d}_start.pdb"
    sysxml = assets / f"{g}_apo_md_rep{rep:02d}_system.xml"
    for p in (start, sysxml):
        if not p.exists():
            raise SystemExit(f"missing {p}")

    pdb = app.PDBFile(str(start))
    system = mm.XmlSerializer.deserialize(sysxml.read_text())
    system.addForce(mm.MonteCarloBarostat(1.0 * u.bar, 300 * u.kelvin, 25))
    integ = mm.LangevinMiddleIntegrator(300 * u.kelvin, 1.0 / u.picosecond,
                                        dt_fs * u.femtoseconds)
    seed = (abs(hash((g, rep))) % 2_000_000_000) + 1
    integ.setRandomNumberSeed(seed)

    try:
        plat = mm.Platform.getPlatformByName("CUDA")
        props = {"Precision": "mixed"}
    except Exception:
        plat, props = mm.Platform.getPlatformByName("CPU"), {}
    sim = app.Simulation(pdb.topology, system, integ, plat, props)
    sim.context.setPositions(pdb.positions)
    sim.context.setVelocitiesToTemperature(300 * u.kelvin, seed)

    keep = np.array([a.index for a in pdb.topology.atoms()
                     if a.residue.name.strip().upper() not in SOLVENT], dtype=np.int32)
    print(f"{g} rep{rep}: {system.getNumParticles()} atoms, keeping {keep.size}, "
          f"platform {sim.context.getPlatform().getName()}, seed {seed}", flush=True)

    n_frames = int(round(ps_total / frame_ps))
    steps = max(1, int(round(frame_ps * 1000.0 / dt_fs)))
    frames = np.empty((n_frames + 1, keep.size, 3), dtype=np.float32)
    # frame 0 is the prepared structure, before any dynamics
    frames[0] = np.asarray(pdb.positions.value_in_unit(u.angstrom))[keep]

    t0 = time.time()
    for i in range(1, n_frames + 1):
        sim.step(steps)
        pos = sim.context.getState(getPositions=True).getPositions(asNumpy=True)
        frames[i] = np.asarray(pos.value_in_unit(u.angstrom))[keep]
        if i % 50 == 0:
            el = time.time() - t0
            print(f"  {i*frame_ps:6.1f} ps  ({el/60:5.1f} min, "
                  f"eta {el/i*(n_frames-i)/60:5.1f} min)", flush=True)

    out.mkdir(parents=True, exist_ok=True)
    np.save(out / f"{g}_rep{rep:02d}_coords.npy", frames)
    np.save(out / f"{g}_rep{rep:02d}_atomindex.npy", keep)
    print(f"wrote {out}/{g}_rep{rep:02d}_coords.npy  {frames.shape}  "
          f"in {(time.time()-t0)/60:.1f} min", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
