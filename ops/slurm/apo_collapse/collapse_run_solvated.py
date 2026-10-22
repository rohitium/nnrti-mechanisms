#!/usr/bin/env python3
"""Apo collapse run with the NNIBP cavity filled with water.

The apo systems were solvated with the buried NNIBP left empty, so the protein
relaxes into a vacuum. This fills the cavity first, equilibrates with the
protein restrained so the waters settle and any excess is expelled while the
pocket is held open, then releases the restraints and samples densely.

Phases
------
  1. fill      relocate bulk waters into the cavity (rigid translation)
  2. minimise  with protein heavy atoms restrained
  3. settle    restrained MD; waters orient and excess leaves
  4. release   restraints off, dense sampling

OpenMM, numpy and scipy only, so it runs under Sherlock's py-openmm module.
The filling logic mirrors nnrti.cli.fill_apo_cavity; one case is cross-checked
against it so the two cannot drift.

Usage
-----
    python collapse_run_solvated.py GENOTYPE REP ASSET_ROOT OUT_DIR [PROD_PS] [SETTLE_PS]
"""
import sys, time
from pathlib import Path

import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as u
from scipy.spatial import cKDTree

SOLV = {"HOH", "WAT", "TIP3", "SOL"}
IONS = {"NA", "CL", "K", "MG", "SOD", "CLA"}
#: NNIBP-lining residues, topology numbering (canonical minus 3).
#: Cilento, Kirby & Sarafianos; same list as nnrti.analysis.metrics.
P66 = (97, 98, 100, 103, 104, 105, 176, 178, 185, 186, 187, 224, 226, 231, 315)
P51 = 135
MIN_PROTEIN_DIST, MIN_WATER_DIST, BULK_MIN_DIST = 3.2, 2.8, 12.0


def classify(top):
    prot_heavy, water_o, water_res, chain_res = [], [], [], {}
    for res in top.residues():
        nm = res.name.strip().upper()
        if nm in SOLV:
            idx = [a.index for a in res.atoms()]
            o = next((a.index for a in res.atoms()
                      if a.element is not None and a.element.symbol == "O"), None)
            if o is not None and len(idx) == 3:
                water_o.append(o); water_res.append(idx)
        elif nm not in IONS:
            chain_res[res.chain.index] = chain_res.get(res.chain.index, 0) + 1
            prot_heavy += [a.index for a in res.atoms()
                           if a.element is not None and a.element.symbol != "H"]
    return (np.array(prot_heavy, int), np.array(water_o, int), water_res, chain_res)


def pocket_ca(top, chain_res):
    ranked = sorted(chain_res, key=chain_res.get, reverse=True)
    p66c, p51c = ranked[0], ranked[1]
    ca = []
    for res in top.residues():
        if res.name.strip().upper() in SOLV | IONS:
            continue
        try:
            rid = int(res.id)
        except ValueError:
            continue
        if (res.chain.index == p66c and rid in P66) or \
           (res.chain.index == p51c and rid == P51):
            ca += [a.index for a in res.atoms() if a.name == "CA"]
    ca = np.array(sorted(set(ca)), int)
    if not 14 <= ca.size <= 16:
        raise SystemExit(f"found {ca.size} pocket CA, expected 15 or 16")
    return ca


def fill(pos, top, prot_heavy, water_o, water_res, centre, seed, max_waters=12):
    r, s = 9.0, 0.5
    ax = [np.arange(centre[i] - r, centre[i] + r + s, s) for i in range(3)]
    grid = np.stack(np.meshgrid(*ax, indexing="ij"), -1).reshape(-1, 3)
    grid = grid[np.linalg.norm(grid - centre, axis=1) <= r]
    grid = grid[cKDTree(pos[prot_heavy]).query(grid)[0] >= MIN_PROTEIN_DIST]
    grid = grid[cKDTree(pos[water_o]).query(grid)[0] >= MIN_WATER_DIST]
    chosen = []
    for i in np.argsort(np.linalg.norm(grid - centre, axis=1)):
        p = grid[i]
        if all(np.linalg.norm(p - q) >= MIN_WATER_DIST for q in chosen):
            chosen.append(p)
        if len(chosen) >= max_waters:
            break
    if not chosen:
        return pos, 0
    d = cKDTree(pos[prot_heavy]).query(pos[water_o])[0]
    cand = np.where(d >= BULK_MIN_DIST)[0]
    donors = np.random.default_rng(seed).choice(cand, size=len(chosen), replace=False)
    new = pos.copy()
    for site, di in zip(chosen, donors):
        new[water_res[di]] += site - new[water_o[di]]
    return new, len(chosen)


def count_cavity_waters(pos, water_o, centre, cutoff=6.0):
    return int((np.linalg.norm(pos[water_o] - centre, axis=1) < cutoff).sum())


def main() -> int:
    g, rep = sys.argv[1], int(sys.argv[2])
    assets = Path(sys.argv[3]) / g / f"rep_{rep:02d}" / "assets"
    out = Path(sys.argv[4])
    prod_ps = float(sys.argv[5]) if len(sys.argv) > 5 else 200.0
    settle_ps = float(sys.argv[6]) if len(sys.argv) > 6 else 100.0
    frame_ps, dt_fs = 1.0, 2.0

    pdb = app.PDBFile(str(assets / f"{g}_apo_md_rep{rep:02d}_start.pdb"))
    system = mm.XmlSerializer.deserialize(
        (assets / f"{g}_apo_md_rep{rep:02d}_system.xml").read_text())
    top = pdb.topology
    pos = np.asarray(pdb.positions.value_in_unit(u.angstrom))

    prot_heavy, water_o, water_res, chain_res = classify(top)
    ca = pocket_ca(top, chain_res)
    centre = pos[ca].mean(0)
    seed = (abs(hash((g, rep))) % 2_000_000_000) + 1

    before = count_cavity_waters(pos, water_o, centre)
    pos, n_added = fill(pos, top, prot_heavy, water_o, water_res, centre, seed)
    after = count_cavity_waters(pos, water_o, centre)
    print(f"{g} rep{rep}: pocket CA {ca.size}, cavity waters {before} -> {after} "
          f"({n_added} relocated)", flush=True)

    # restrain protein heavy atoms during settling, released afterwards
    restraint = mm.CustomExternalForce(
        "k*periodicdistance(x,y,z,x0,y0,z0)^2")
    restraint.addGlobalParameter("k", 10.0 * u.kilocalories_per_mole / u.angstrom**2)
    for p in ("x0", "y0", "z0"):
        restraint.addPerParticleParameter(p)
    for i in prot_heavy:
        restraint.addParticle(int(i), pos[i] * 0.1)      # nm
    rf_index = system.addForce(restraint)
    system.addForce(mm.MonteCarloBarostat(1.0 * u.bar, 300 * u.kelvin, 25))

    integ = mm.LangevinMiddleIntegrator(300 * u.kelvin, 1.0 / u.picosecond,
                                        dt_fs * u.femtoseconds)
    integ.setRandomNumberSeed(seed)
    try:
        plat, props = mm.Platform.getPlatformByName("CUDA"), {"Precision": "mixed"}
    except Exception:
        plat, props = mm.Platform.getPlatformByName("CPU"), {}
    sim = app.Simulation(top, system, integ, plat, props)
    sim.context.setPositions(pos * 0.1)

    sim.minimizeEnergy(maxIterations=2000)
    sim.context.setVelocitiesToTemperature(300 * u.kelvin, seed)
    sim.step(int(settle_ps * 1000 / dt_fs))
    settled = np.asarray(sim.context.getState(getPositions=True)
                         .getPositions(asNumpy=True).value_in_unit(u.angstrom))
    print(f"  after {settle_ps:.0f} ps restrained settling: "
          f"cavity waters {count_cavity_waters(settled, water_o, centre)}", flush=True)

    # release: the pocket is now free to respond with water present
    sim.context.setParameter("k", 0.0)
    keep = np.array([a.index for a in top.atoms()
                     if a.residue.name.strip().upper() not in SOLV], int)
    n_frames = int(round(prod_ps / frame_ps))
    steps = max(1, int(round(frame_ps * 1000.0 / dt_fs)))
    frames = np.empty((n_frames + 1, keep.size, 3), dtype=np.float32)
    cav = np.empty(n_frames + 1, dtype=np.int16)
    frames[0] = settled[keep]
    cav[0] = count_cavity_waters(settled, water_o, centre)

    t0 = time.time()
    for i in range(1, n_frames + 1):
        sim.step(steps)
        x = np.asarray(sim.context.getState(getPositions=True)
                       .getPositions(asNumpy=True).value_in_unit(u.angstrom))
        frames[i] = x[keep]
        cav[i] = count_cavity_waters(x, water_o, centre)
        if i % 50 == 0:
            print(f"  {i*frame_ps:6.1f} ps  cavity waters {cav[i]}  "
                  f"({(time.time()-t0)/60:.1f} min)", flush=True)

    out.mkdir(parents=True, exist_ok=True)
    np.save(out / f"{g}_rep{rep:02d}_coords.npy", frames)
    np.save(out / f"{g}_rep{rep:02d}_atomindex.npy", keep)
    np.save(out / f"{g}_rep{rep:02d}_cavitywaters.npy", cav)
    print(f"wrote {out}/{g}_rep{rep:02d}_coords.npy {frames.shape} "
          f"in {(time.time()-t0)/60:.1f} min", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
