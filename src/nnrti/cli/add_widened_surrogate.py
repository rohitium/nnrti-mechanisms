"""Fill in surrogate columns that the candidate scan does not produce.

stage_scan measures only the partners in surrogate_candidates.csv, which the
cylinder filter (PERP_MAX / REACH_MIN / REACH_MAX) keeps close to the
protein-ligand contact direction. The residue-106 surrogate came out of a
deliberately widened search instead -- every backbone partner 6-20 A from the
residue, 399 of them against the 9 the cylinder allows -- and so its column is
absent whenever the scan is re-run. That used to mean re-deriving it by hand,
which is how it was lost once already.

This reads the chosen surrogates, works out which of their columns the
per-frame table is missing, and measures those with the same anchors, the same
minimum-image convention and the same trailing window as the scan.
"""

from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd

from nnrti.cli.compute_mechanism_coordinates import OFFSET, _replicate_inputs
from nnrti.cli.find_pocket_surrogates import (
    AXES,
    REPO_ROOT,
    _anchor,
    _apo_replicates,
    _group_centroid,
    _mic,
    _res_atoms,
    _run_ns,
    _window_mask,
)

KEYS = ["mutation", "state", "replicate", "frame"]

# "d(106...DOR)|sc|ASP186:CA" -> axis, anchor tag, partner residue and atom
COLUMN = re.compile(r"^(?P<axis>.+?)\|(?P<tag>[^|]+)\|(?P<resname>[A-Z0-9]+)"
                    r"(?P<canonical>\d+):(?P<atom>\S+)$")


def _anchor_group(top, tag: str, canonical: int, names):
    """The atom indices a surrogate's anchor tag refers to.

    "sc" is the side chain alone, which is what the Table 3 axis for residues
    103 and 106 measures from; "ca" is the single CA, which cannot move because
    a side chain was substituted; "as_table3" is whatever atom set the axis
    itself uses.
    """
    if tag == "ca":
        return _res_atoms(top, canonical, names={"CA"})
    if tag == "sc":
        return _anchor(top, canonical, None)
    if tag == "as_table3":
        return _anchor(top, canonical, names)
    raise ValueError(f"unknown anchor tag: {tag}")


def main(argv: list[str] | None = None) -> int:
    import mdtraj as md

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest", type=Path, default=Path("manifests/md_manifest.csv"))
    p.add_argument("--output-dir", type=Path,
                   default=Path("results/analysis/pocket_surrogates"))
    p.add_argument("--window-ns", type=float, default=50.0)
    p.add_argument("--min-ns", type=float, default=80.0)
    args = p.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    tables = args.output_dir / "tables"
    per_frame_path = tables / "surrogate_per_frame_both.csv"
    frames = pd.read_csv(per_frame_path)
    selected = pd.read_csv(tables / "surrogates_selected.csv")

    wanted = [s for s in selected["surrogate"]
              if s != "DIRECT" and s not in frames.columns]
    if not wanted:
        logging.info("every selected surrogate is already measured; nothing to add")
        return 0
    specs = {axis: names for axis, _c, _a, names, _lg in AXES}
    anchors = {axis: canonical for axis, _c, canonical, _n, _lg in AXES}

    parsed = []
    for column in wanted:
        m = COLUMN.match(column)
        if m is None:
            logging.error("cannot parse surrogate column: %s", column)
            return 1
        axis = m.group("axis")
        if axis not in anchors:
            logging.error("%s is not one of the Table 3 axes", axis)
            return 1
        parsed.append(dict(column=column, axis=axis, tag=m.group("tag"),
                           canonical_anchor=anchors[axis],
                           partner_resid=int(m.group("canonical")) + OFFSET,
                           partner_atom=m.group("atom")))
        logging.info("adding %s", column)

    mf = pd.read_csv(args.manifest)
    jobs = []
    for _, r in mf.iterrows():
        try:
            topo, dcd, _ = _replicate_inputs(r, REPO_ROOT)
            jobs.append(("holo", str(r["mutation"]), int(r["replicate"]), topo, dcd))
        except Exception as exc:
            logging.error("holo %s rep%s: %s", r["mutation"], r["replicate"], exc)
    for mutation in dict.fromkeys(str(m) for m in mf["mutation"]):
        for rep, topo, dcd in _apo_replicates(mutation):
            jobs.append(("apo", mutation, rep, topo, dcd))

    rows = []
    for state, mutation, rep, topo, dcd in jobs:
        try:
            traj = md.load(str(dcd), top=str(topo))
            total_ns = _run_ns(Path(dcd).parent)
            if total_ns < args.min_ns:
                continue
            if traj.unitcell_lengths is None:
                raise ValueError("no box information; minimum image is undefined")
            top = traj.topology
            xyz = traj.xyz * 10.0
            box = traj.unitcell_lengths * 10.0

            series = {}
            for spec in parsed:
                idx = _anchor_group(top, spec["tag"], spec["canonical_anchor"],
                                    specs[spec["axis"]])
                partner = [a.index for a in top.atoms
                           if a.residue.resSeq == spec["partner_resid"]
                           and a.name == spec["partner_atom"]]
                if not idx or not partner:
                    raise ValueError(f"no atoms for {spec['column']}")
                centre = _group_centroid(xyz, list(idx), box)
                series[spec["column"]] = np.linalg.norm(
                    _mic(centre - xyz[:, partner[0], :], box), axis=1)

            mask = _window_mask(traj.n_frames, total_ns, args.window_ns)
            for f in np.flatnonzero(mask):
                f = int(f)
                rec = dict(mutation=mutation, state=state, replicate=rep, frame=f)
                for column, values in series.items():
                    rec[column] = float(values[f])
                rows.append(rec)
            logging.info("%s %s rep%02d: %d frames", mutation, state, rep,
                         int(mask.sum()))
        except Exception as exc:
            logging.error("FAILED %s %s rep%s: %s", mutation, state, rep, exc)

    if not rows:
        logging.error("no frames measured")
        return 1
    added = pd.DataFrame(rows)
    merged = frames.merge(added, on=KEYS, how="left", validate="one_to_one")
    for column in wanted:
        filled = int(merged[column].notna().sum())
        logging.info("%s: %d of %d rows populated", column, filled, len(merged))
        if filled != len(merged):
            logging.error("%s is incomplete; refusing to overwrite the table", column)
            return 1
    merged.to_csv(per_frame_path, index=False)
    logging.info("wrote %s", per_frame_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
