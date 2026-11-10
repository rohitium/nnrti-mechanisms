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

#: The crystal structure uses canonical HIV-1 RT numbering, so no offset. The
#: simulation topology is shifted by -3, which --resid-offset handles if the
#: figure is ever rebuilt from an MD frame instead.
RESID_OFFSET = 0

#: NNIBP-lining residues, canonical numbering (Cilento, Kirby & Sarafianos).
NNIBP_P66 = (100, 101, 103, 106, 107, 108, 179, 181, 188, 189, 190, 227, 229, 234, 318)
NNIBP_P51 = (138,)

#: Chain IDs. In 4NCG the ligand sits in chain A with p66, so DOR is selected
#: by residue name rather than by chain.
CHAIN_P66, CHAIN_P51 = "A", "B"
LIGAND_RESN = "2KW"


PYMOL_APP = Path("/Applications/PyMOL.app/Contents/MacOS/PyMOL")


def mutated_positions(manifest: Path, offset: int) -> pd.DataFrame:
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
             "topology": p + offset,
             "substitutions": ",".join(sorted(e["subs"])),
             "max_fold_reduction": e["max_fold"],
             "in_nnibp": p in NNIBP_P66}
            for p, e in sorted(seen.items())]
    return pd.DataFrame(rows)




def prepare_structure(src: Path, out_dir: Path, pymol: Path) -> Path:
    """Return a solvent-free PDB of `src`, converting from mmCIF if needed.

    Everything downstream (mdtraj for the view search, PyMOL for the scene)
    then reads one file in one format, instead of each needing its own mmCIF
    support.
    """
    if src.suffix.lower() not in (".cif", ".mmcif"):
        return src
    dest = out_dir / f"{src.stem}_nosolvent.pdb"
    script = out_dir / "_convert.pml"
    script.write_text(
        f"load {src.resolve()}, s\n"
        f"remove s and solvent\n"
        f"remove s and not (polymer or resn {LIGAND_RESN})\n"
        f"save {dest.resolve()}, s\n")
    r = subprocess.run([str(pymol), "-cq", str(script)],
                       capture_output=True, text=True)
    script.unlink(missing_ok=True)
    if r.returncode != 0 or not dest.exists():
        print(r.stdout[-1500:], r.stderr[-1500:])
        raise SystemExit(f"could not convert {src}")
    return dest


def detect_offset(structure: Path) -> int:
    """Residue numbering offset from canonical, read off the structure.

    Canonical position 225 is a proline in HIV-1 RT, so whichever numbering
    puts PRO at 225 is the canonical one. Guessing wrong would silently select
    the wrong residues, so this checks rather than assumes.
    """
    import mdtraj as md
    top = md.load(str(structure)).topology
    for offset in (0, -3):
        hit = [r for r in top.chain(0).residues if r.resSeq == 225 + offset]
        if hit and hit[0].name == "PRO":
            return offset
    raise SystemExit(f"{structure}: cannot find PRO225 at offset 0 or -3; "
                     "pass --resid-offset explicitly")


def best_view(structure: Path, muts: pd.DataFrame,
              n_directions: int = 4000) -> tuple[np.ndarray, np.ndarray, float]:
    """Camera rotation that separates the labelled residues best on screen.

    Labels collide when two mutated CA atoms project close together, which is a
    property of the viewing direction and can be optimised rather than hunted
    for by hand. Directions are sampled on a Fibonacci sphere; for each, the
    twelve CA atoms are projected onto the image plane and scored by the
    smallest pairwise separation. A direction that puts a label on top of the
    ligand is penalised, since the drug has to stay readable too.

    Returns the rotation matrix (world to camera, row-major), the centre to
    rotate about, and the radius the camera needs to clear.
    """
    import mdtraj as md

    t = md.load(str(structure))
    top = t.topology
    want = {r.topology for r in muts.itertuples()}
    ca = [a.index for a in top.atoms
          if a.name == "CA" and a.residue.chain.index == 0
          and a.residue.resSeq in want]
    lig = [a.index for a in top.atoms if a.residue.name == "2KW"]
    if len(ca) != len(muts) or not lig:
        raise SystemExit(f"expected {len(muts)} mutated CA atoms and a ligand; "
                         f"found {len(ca)} and {len(lig)}")

    xyz = t.xyz[0] * 10.0
    pts = xyz[ca]
    lig_c = xyz[lig].mean(axis=0)
    centre = np.vstack([pts, xyz[lig]]).mean(axis=0)
    radius = float(np.linalg.norm(np.vstack([pts, xyz[lig]]) - centre, axis=1).max())

    # Fibonacci sphere: even coverage of directions without a pole bias.
    i = np.arange(n_directions) + 0.5
    phi = np.arccos(1.0 - 2.0 * i / n_directions)
    theta = np.pi * (1.0 + 5.0 ** 0.5) * i
    dirs = np.c_[np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi),
                 np.cos(phi)]

    best, best_R = -np.inf, None
    for d in dirs:
        z = d / np.linalg.norm(d)
        up = np.array([0.0, 0.0, 1.0])
        if abs(np.dot(up, z)) > 0.9:
            up = np.array([0.0, 1.0, 0.0])
        x = np.cross(up, z); x /= np.linalg.norm(x)
        y = np.cross(z, x)
        R = np.vstack([x, y, z])

        proj = (pts - centre) @ R[:2].T
        dists = np.linalg.norm(proj[:, None, :] - proj[None, :, :], axis=-1)
        np.fill_diagonal(dists, np.inf)
        min_sep = float(dists.min())

        # Keep labels off the ligand: penalise any that lands within 4 A of it
        # in projection, so the drug stays legible under its own labels.
        lig_proj = (lig_c - centre) @ R[:2].T
        clearance = float(np.linalg.norm(proj - lig_proj, axis=1).min())
        score = min_sep + min(clearance, 6.0)

        if score > best:
            best, best_R = score, R
    print(f"view search: {n_directions} directions, best label separation "
          f"score {best:.1f} A")
    return best_R, centre, radius, pts


def write_pml(pml: Path, structure: Path, muts: pd.DataFrame,
              offset: int, view=None) -> None:
    nn66 = "+".join(str(r + offset) for r in NNIBP_P66)
    nn51 = "+".join(str(r + offset) for r in NNIBP_P51)
    mut_resi = "+".join(str(r.topology) for r in muts.itertuples())

    L: list[str] = [
        "# Figure 1: mutations mapped onto RT-DOR.",
        "# Built by nnrti.cli.build_figure1_pymol -- re-run it to regenerate.",
        "# Everything below is a starting point; adjust and re-save the session.",
        "",
        "reinitialize",
        f"load {structure}, RT_DOR",
        # Crystallographic waters are not part of this figure.
        "remove RT_DOR and solvent",
        "hide everything",
        "",
        "# ---- the three molecules -------------------------------------------",
        f"select p66, RT_DOR and chain {CHAIN_P66} and polymer",
        f"select p51, RT_DOR and chain {CHAIN_P51} and polymer",
        f"select DOR, RT_DOR and resn {LIGAND_RESN}",
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
          "color grey70, p51",
          # The surrounding fold is context, not subject. Fading it well back
          # leaves DOR and the mutated side chains carrying the image, as in
          # the TOC figure.
          "set cartoon_transparency, 0.82, p66",
          "set cartoon_transparency, 0.90, p51",
          "set cartoon_transparency, 0.45, NNIBP",
          "",
          # DOR: yellow carbons with explicit CPK heteroatoms. Set element by
          # element rather than calling util.cnc, so fluorine is the aqua asked
          # for instead of PyMOL's default pale green and nothing depends on
          # which element palette a given PyMOL build ships.
          "hide everything, DOR and hydro",
          "show sticks, DOR",
          "color yellow,  DOR and elem C",
          "color blue,    DOR and elem N",
          "color red,     DOR and elem O",
          "color aquamarine, DOR and elem F",
          "color green,   DOR and elem Cl",
          "show spheres, DOR",
          "set sphere_scale, 0.20, DOR",
          "",
          # Keep the backbone N for proline only. P225's side chain is a ring
          # that closes through its own N (N-CA-CB-CG-CD-N), so excluding N
          # breaks two bonds and leaves CB-CG-CD dangling off CA. Every other
          # residue drops N, which would otherwise draw a stub into the
          # backbone. Glycine (G190) has no side chain at all and is shown by
          # its CA sphere alone, which is correct rather than missing.
          "show sticks, mutations and not hydro and not (name C+O) "
          "and not (name N and not resn PRO)",
          "show spheres, mutations and name CA",
          # Small enough not to swallow the side chain it marks. At 0.45 the
          # sphere covered proline's CA vertex and the two ring bonds meeting
          # there, so P225 read as an open chain even with every ring atom
          # drawn.
          "set sphere_scale, 0.26, mutations and name CA",
          "",
          "show cartoon, NNIBP",
          "color palecyan, NNIBP",
          "",
          "# ---- colour the mutated sites --------------------------------------",
          "color orange, mutations",
          "color orange, mutations and elem C",
          "color blue,   mutations and elem N",
          "color red,    mutations and elem O",
          "color yellow, mutations and elem S",
    ]

    # The b-factor still carries log10(fold reduction) at each mutated position,
    # so the flat orange above can be swapped for a phenotype scale by
    # uncommenting one line. Nothing else in the session depends on it.
    L += [
        "",
        "# b-factor at each mutated position is log10(fold reduction in DOR",
        "# susceptibility), left in place so the sites can be coloured by",
        "# phenotype instead of flat orange. Uncomment the spectrum line.",
        "alter RT_DOR, b=0.0",
    ]
    for r in muts.itertuples():
        b = float(np.log10(r.max_fold_reduction)) if r.max_fold_reduction > 0 else 0.0
        L.append(f'alter mut_{r.wt_aa}{r.canonical}, b={b:.3f}')
    hi = float(np.log10(muts.max_fold_reduction.max()))
    L += [
        "rebuild",
        f"# spectrum b, yellow_orange_red, mutations, minimum=0.00, maximum={hi:.2f}",
    ]

    L += [
        "# ---- labels ----------------------------------------------------------",
        "# Canonical numbering, which is what the manuscript uses.",
    ]
    for r in muts.itertuples():
        L.append(f'label mut_{r.wt_aa}{r.canonical} and name CA, "{r.wt_aa}{r.canonical}"')
    L += [
        "set label_size, 20",
        # Black on the white background, with a white halo so a label crossing
        # a dark side chain stays readable.
        "set label_color, black",
        "set label_outline_color, white",
        "set label_shadow_mode, 0",
        "set label_bg_color, white",
        "set label_bg_transparency, 0.25",
        "set label_position, (0, 0, 3.0)",
        "",
        "# ---- label nudges: edit these ----------------------------------------",
        "# One line per residue. The three numbers are a screen-space offset in",
        "# angstroms: x right, y up, z toward the viewer. The anchor sphere does",
        "# not move, so a nudged label still marks its own residue.",
        "#",
        "# In the GUI you can also drag a label directly: enable editing with",
        "# the Mouse menu (3-Button Editing), then ctrl-shift-drag the label.",
        "# PyMOL stores the result in these same settings, so File > Save",
        "# Session keeps whatever you do by hand.",
        "",
        "# ---- scene -----------------------------------------------------------",
        "bg_color white",
        "set ray_opaque_background, 1",
        "set cartoon_fancy_helices, 1",
        "set antialias, 2",
        "set ray_shadows, 0",
        "deselect",
        "",
        "# Suggested extras, left off so the first view stays readable:",
        "#   show surface, NNIBP",
        "#   set transparency, 0.5",
        "#   show cartoon, p51   (p51 is already shown, semi-transparent)",
    ]
    for r in muts.itertuples():
        L.append(f'set label_position, (0, 0, 3.0), mut_{r.wt_aa}{r.canonical}')

    if view is not None:
        R, centre, radius, _ = view
        # PyMOL's default field of view is 20 degrees, so a sphere of this
        # radius only fits the frame beyond radius / tan(10 deg) = 5.7 x. The
        # first attempt used 3.2 x and cropped A98, M230 and Y318 off-screen.
        dist = 6.6 * radius
        m = [f"{v:.6f}" for v in R.flatten()]
        o = [f"{v:.6f}" for v in centre]
        L += [
            "",
            "# Camera chosen by best_view(): of 4000 sampled directions, the",
            "# one that maximises the smallest on-screen separation between",
            "# labelled residues while keeping them clear of the ligand. It is",
            "# a starting point, not a final answer -- P225/F227 and A98/Y318",
            "# are only a few angstroms apart, so no direction separates them.",
            "# Rotate as you like and re-save; nothing else depends on it.",
            "set_view (\\",
            f"  {m[0]}, {m[1]}, {m[2]},\\",
            f"  {m[3]}, {m[4]}, {m[5]},\\",
            f"  {m[6]}, {m[7]}, {m[8]},\\",
            f"  0.0, 0.0, {-dist:.3f},\\",
            f"  {o[0]}, {o[1]}, {o[2]},\\",
            f"  {dist - radius:.3f}, {dist + radius:.3f}, -20.0 )",
        ]
    else:
        L += ["orient (DOR or mutations)", "zoom (DOR or mutations), 6"]
    pml.write_text("\n".join(L) + "\n")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--structure", type=Path,
                   default=Path("data/structures/4NCG.cif"),
                   help="WT RT-DOR structure. Defaults to the 4NCG crystal "
                        "structure, which is what Figure 1 shows.")
    p.add_argument("--resid-offset", type=int, default=None,
                   help="Residue numbering offset from canonical. 0 for a "
                        "crystal structure, -3 for the simulation topology. "
                        "Detected from the structure when not given.")
    p.add_argument("--manifest", type=Path, default=Path("manifests/md_manifest.csv"))
    p.add_argument("--output-dir", type=Path, default=Path("paper/generated/figure1"))
    p.add_argument("--pymol", type=Path, default=PYMOL_APP)
    args = p.parse_args(argv)

    if not args.structure.exists():
        raise SystemExit(f"structure not found: {args.structure}")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    structure = prepare_structure(args.structure, args.output_dir, args.pymol)
    offset = args.resid_offset if args.resid_offset is not None \
        else detect_offset(structure)
    print(f"{args.structure.name}: residue numbering offset {offset:+d} "
          f"({'canonical' if offset == 0 else 'simulation topology'})")
    muts = mutated_positions(args.manifest, offset)
    pml = args.output_dir / "figure1_mutations.pml"
    pse = args.output_dir / "figure1_mutations.pse"
    table = args.output_dir / "figure1_mutated_positions.csv"
    muts.to_csv(table, index=False)
    view = best_view(structure, muts)
    write_pml(pml, structure.resolve(), muts, offset, view)

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
