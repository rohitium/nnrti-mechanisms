#!/usr/bin/env python3
"""Figure 1: every studied mutation mapped onto the RT-DOR structure.

Reviewer 1 (comment 4) asked for a figure showing where the mutations sit
relative to DOR and the NNIBP, because the structural information is otherwise
spread across dense tables.

This writes a PyMOL script and runs it to produce a .pse session rather than a
finished image. Orientation, label placement and ray-tracing are judgement
calls best made by eye, so the session is built to be opened and adjusted:
every group of atoms is a named selection, nothing is flattened into a single
representation, and the view is only a starting point.

Residue numbering: the simulation topology is offset from the canonical HIV-1
RT numbering the manuscript uses by RESID_OFFSET. Selections are written in
topology numbering so they resolve, and every label shows the canonical number,
so what the reader sees matches the text.

Mutated positions and their phenotypes are read from the manifest rather than
hard-coded, so the figure cannot drift from the panel actually simulated.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.build_figure1_pymol
"""

from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

#: Simulation topology numbering = canonical numbering + RESID_OFFSET.
RESID_OFFSET = -3

#: NNIBP-lining residues, canonical numbering (Cilento, Kirby & Sarafianos).
NNIBP_P66 = (100, 101, 103, 106, 107, 108, 179, 181, 188, 189, 190, 227, 229, 234, 318)
NNIBP_P51 = (138,)

#: Chain IDs in the analysis topology.
CHAIN_P66, CHAIN_P51, CHAIN_DOR = "A", "B", "E"

PYMOL_APP = Path("/Applications/PyMOL.app/Contents/MacOS/PyMOL")


def mutated_positions(manifest: Path) -> pd.DataFrame:
    """One row per mutated position, with the worst phenotype seen there."""
    m = pd.read_csv(manifest)
    seen: dict[int, dict] = {}
    for _, r in m[["mutation", "fold_reduction"]].drop_duplicates().iterrows():
        if str(r.mutation).upper() == "WT":
            continue
        fold = float(r.fold_reduction) if pd.notna(r.fold_reduction) else np.nan
        for token in str(r.mutation).split("+"):
            hit = re.match(r"([A-Z])(\d+)([A-Z])", token.strip())
            if not hit:
                continue
            wt_aa, pos = hit.group(1), int(hit.group(2))
            e = seen.setdefault(pos, {"wt_aa": wt_aa, "subs": set(), "max_fold": np.nan})
            e["subs"].add(token.strip())
            e["max_fold"] = np.nanmax([e["max_fold"], fold])
    rows = [{"canonical": p, "wt_aa": e["wt_aa"],
             "topology": p + RESID_OFFSET,
             "substitutions": ",".join(sorted(e["subs"])),
             "max_fold_reduction": e["max_fold"],
             "in_nnibp": p in NNIBP_P66}
            for p, e in sorted(seen.items())]
    return pd.DataFrame(rows)


def write_pml(pml: Path, structure: Path, muts: pd.DataFrame) -> None:
    nn66 = "+".join(str(r + RESID_OFFSET) for r in NNIBP_P66)
    nn51 = "+".join(str(r + RESID_OFFSET) for r in NNIBP_P51)
    mut_resi = "+".join(str(r.topology) for r in muts.itertuples())

    L: list[str] = [
        "# Figure 1: mutations mapped onto RT-DOR.",
        "# Built by nnrti.cli.build_figure1_pymol -- re-run it to regenerate.",
        "# Everything below is a starting point; adjust and re-save the session.",
        "",
        "reinitialize",
        f'load {structure}, RT_DOR',
        "hide everything",
        "",
        "# ---- the three molecules -------------------------------------------",
        f'select p66, RT_DOR and chain {CHAIN_P66} and polymer',
        f'select p51, RT_DOR and chain {CHAIN_P51} and polymer',
        f'select DOR, RT_DOR and chain {CHAIN_DOR}',
        "",
        "# ---- the pocket ------------------------------------------------------",
        f'select NNIBP, (p66 and resi {nn66}) or (p51 and resi {nn51})',
        "",
        "# ---- mutated positions ----------------------------------------------",
        f'select mutations, p66 and resi {mut_resi}',
    ]

    # One selection per position, named with the canonical number the paper uses.
    for r in muts.itertuples():
        L.append(f'select mut_{r.wt_aa}{r.canonical}, p66 and resi {r.topology}')
    L += ["", "# ---- representations -------------------------------------------------",
          "show cartoon, p66 or p51",
          "color grey80, p66",
          "color grey60, p51",
          "set cartoon_transparency, 0.55, p51",
          "",
          "show sticks, DOR",
          "color yellow, DOR and elem C",
          "util.cnc('DOR')",
          "show spheres, DOR",
          "set sphere_scale, 0.22, DOR",
          "",
          "show sticks, mutations and not (name C+N+O)",
          "show spheres, mutations and name CA",
          "set sphere_scale, 0.45, mutations and name CA",
          "",
          "show cartoon, NNIBP",
          "color palecyan, NNIBP",
          "",
          "# ---- colour the mutated sites by phenotype severity -----------------",
          "# b-factor is overwritten with log10(max fold reduction) at each mutated",
          "# position, so the spectrum below is a phenotype scale, not a thermal one.",
          "alter RT_DOR, b=0.0",
    ]
    for r in muts.itertuples():
        b = float(np.log10(r.max_fold_reduction)) if r.max_fold_reduction > 0 else 0.0
        L.append(f'alter mut_{r.wt_aa}{r.canonical}, b={b:.3f}')
    lo, hi = 0.0, float(np.log10(muts.max_fold_reduction.max()))
    L += [
        "rebuild",
        f"spectrum b, white_yellow_orange_red, mutations, minimum={lo:.2f}, maximum={hi:.2f}",
        "",
        "# ---- labels ----------------------------------------------------------",
        "# Canonical numbering, which is what the manuscript uses.",
    ]
    for r in muts.itertuples():
        L.append(f'label mut_{r.wt_aa}{r.canonical} and name CA, "{r.wt_aa}{r.canonical}"')
    L += [
        "set label_size, 16",
        "set label_color, black",
        "set label_outline_color, white",
        "set label_position, (0, 0, 2.5)",
        "",
        "# ---- scene -----------------------------------------------------------",
        "bg_color white",
        "set ray_opaque_background, 1",
        "set cartoon_fancy_helices, 1",
        "set antialias, 2",
        "set ray_shadows, 0",
        "orient (DOR or mutations)",
        "zoom (DOR or mutations), 6",
        "deselect",
        "",
        "# Suggested extras, left off so the first view stays readable:",
        "#   show surface, NNIBP",
        "#   set transparency, 0.5",
        "#   show cartoon, p51   (p51 is already shown, semi-transparent)",
    ]
    pml.write_text("\n".join(L) + "\n")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--structure", type=Path,
                   default=Path("results/md_runs/wt/rep_01/wt_rep01_analysis_topology.pdb"),
                   help="WT RT-DOR structure. Defaults to the analysis topology, "
                        "i.e. the system actually simulated, stripped of solvent.")
    p.add_argument("--manifest", type=Path, default=Path("manifests/md_manifest.csv"))
    p.add_argument("--output-dir", type=Path, default=Path("revisions/figure1"))
    p.add_argument("--pymol", type=Path, default=PYMOL_APP)
    args = p.parse_args(argv)

    if not args.structure.exists():
        raise SystemExit(f"structure not found: {args.structure}")
    muts = mutated_positions(args.manifest)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    pml = args.output_dir / "figure1_mutations.pml"
    pse = args.output_dir / "figure1_mutations.pse"
    table = args.output_dir / "figure1_mutated_positions.csv"
    muts.to_csv(table, index=False)
    write_pml(pml, args.structure.resolve(), muts)

    print(f"{len(muts)} mutated positions, "
          f"{int(muts.in_nnibp.sum())} of them NNIBP-lining:")
    print(muts.to_string(index=False))

    if not args.pymol.exists():
        print(f"\nPyMOL not found at {args.pymol}; wrote {pml} only. "
              f"Open it with:  pymol {pml}")
        return 0

    script = f"{pml.resolve()}\nsave {pse.resolve()}\n"
    runner = args.output_dir / "_build_session.pml"
    runner.write_text(f"@{pml.resolve()}\nsave {pse.resolve()}\n")
    r = subprocess.run([str(args.pymol), "-cq", str(runner)],
                       capture_output=True, text=True)
    runner.unlink(missing_ok=True)
    if r.returncode != 0 or not pse.exists():
        print(r.stdout[-2000:])
        print(r.stderr[-2000:])
        raise SystemExit("PyMOL failed to build the session")
    print(f"\nwrote {pml}\n      {pse}\n      {table}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
