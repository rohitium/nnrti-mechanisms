#!/usr/bin/env python3
"""Table-of-contents graphic: render in PyMOL, compose and label at print size.

ACS wants the TOC graphic on the manuscript's last page as .tif or .pdf, with
no caption and all text legible at roughly 3.25 x 1.75 inches.

Two stages, deliberately separated. PyMOL ray-traces the molecular scene
(ops/figures/toc_graphic.pml); this module then places that render on a canvas
of the exact printed size and adds the text in points. Baking labels into the
ray trace ties their size to whatever pixel dimensions the render happened to
use, which is the usual reason TOC text comes out unreadable.

The figure shows one thing -- doravirine bound in the NNIBP, with the panel's
substituted positions around it -- and states the result in words. An earlier
version paired a projected Ca trace with the ddG-against-phenotype scatter.
The trace read as scribble rather than structure, and leading with the
scatter put the weakest result in the paper in front of a reader deciding
whether to open it.

Usage
-----
    PYTHONPATH=src python -m nnrti.cli.plot_toc_graphic
    PYTHONPATH=src python -m nnrti.cli.plot_toc_graphic --skip-render
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

PYMOL = Path("/Applications/PyMOL.app/Contents/bin/pymol")
WIDTH_IN, HEIGHT_IN = 3.3, 1.8


def trim_white(path: Path, pad: int = 8, drop_title: bool = True) -> Path:
    """Trim the schematic's white border, and drop its title band.

    Two reasons. The source render carries generous margins, which would shrink
    the cycle to the middle of its panel and waste the space. And the figure's
    own "Thermodynamic cycle" heading is redundant here: a JCIM reader
    recognises the construction on sight, and the TOC graphic takes no caption,
    so the words only cost space the diagram can use.

    The title is found structurally rather than by pixel offset: rows are
    grouped into bands of ink separated by blank gaps, and a short leading band
    is the heading. That survives a re-render at a different size.
    """
    import numpy as np
    from PIL import Image

    im = Image.open(path).convert("RGB")
    a = np.asarray(im)
    ink = (a < 245).any(axis=2)
    rows, cols = np.where(ink.any(axis=1))[0], np.where(ink.any(axis=0))[0]
    if not len(rows) or not len(cols):
        return path

    top = int(rows[0])
    if drop_title:
        filled = ink.any(axis=1)
        bands, run = [], None
        for y, on in enumerate(filled):
            if on and run is None:
                run = y
            elif not on and run is not None:
                bands.append((run, y - 1))
                run = None
        if run is not None:
            bands.append((run, len(filled) - 1))
        # A short first band with a real gap after it is the heading.
        if len(bands) >= 2:
            h0 = bands[0][1] - bands[0][0]
            gap = bands[1][0] - bands[0][1]
            if h0 < 0.15 * im.height and gap > 0.02 * im.height:
                top = bands[1][0]

    box = (max(int(cols[0]) - pad, 0), max(top - pad, 0),
           min(int(cols[-1]) + pad, im.width), min(int(rows[-1]) + pad, im.height))
    out = Path("/tmp/toc_cycle_trimmed.png")
    im.crop(box).save(out)
    return out


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pml", type=Path, default=Path("ops/figures/toc_graphic.pml"))
    p.add_argument("--scene", type=Path, default=Path("/tmp/toc_scene.png"))
    p.add_argument("--cycle", type=Path,
                   default=Path("results/analysis/fep_pmx/protocol_schematic/"
                                "01_thermodynamic_cycle.png"),
                   help="Thermodynamic cycle for the right-hand panel. The original "
                        "200 dpi render, not a crop out of the figure PDF, which "
                        "rasterises badly at TOC size.")
    p.add_argument("--output-dir", type=Path, default=Path("revisions"))
    p.add_argument("--dpi", type=int, default=600)
    p.add_argument("--skip-render", action="store_true",
                   help="Reuse the existing scene render instead of invoking PyMOL.")
    args = p.parse_args(argv)

    if not args.skip_render:
        if not PYMOL.exists():
            raise SystemExit(f"PyMOL not found at {PYMOL}; pass --skip-render to reuse "
                             f"{args.scene}")
        subprocess.run([str(PYMOL), "-cq", str(args.pml)], check=True)
    if not args.scene.exists():
        raise SystemExit(f"no scene render at {args.scene}")

    img = mpimg.imread(args.scene)

    fig = plt.figure(figsize=(WIDTH_IN, HEIGHT_IN))
    fig.patch.set_facecolor("white")

    # Left: the molecular scene. Right: Figure 2 of the paper, the alchemical
    # free energy protocol. The protocol is the point of the work -- it is what
    # makes the binding free energy an equilibrium quantity rather than an
    # end-point score -- so it earns the second half of the graphic. MM/GBSA is
    # a secondary comparison in the paper and does not belong here.
    ax = fig.add_axes([0.005, 0.02, 0.43, 0.96])
    ax.imshow(img, interpolation="lanczos")
    ax.axis("off")

    fig2 = mpimg.imread(str(trim_white(args.cycle)))
    h, w = fig2.shape[:2]
    box_w = 0.555
    box_h = box_w * (h / w) * (WIDTH_IN / HEIGHT_IN)
    ax2 = fig.add_axes([0.44, 0.575 - box_h / 2, box_w, box_h])
    ax2.imshow(fig2, interpolation="lanczos")
    ax2.axis("off")

    # The one label on the graphic. Bottom left over the structure, on a single
    # line: "mutated" rather than "substituted" so it matches the Mut RT states
    # in the cycle beside it.
    ax.scatter([0.045], [0.055], s=30, color="#f08c1b", edgecolor="#6b3300",
               linewidths=0.45, transform=ax.transAxes, zorder=9, clip_on=False)
    ax.text(0.095, 0.055, "mutated positions", transform=ax.transAxes,
            fontsize=6.0, color="#2e2e2e", va="center", ha="left", zorder=9)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(args.output_dir / f"TOC-Graphic.{ext}", dpi=args.dpi,
                    facecolor="white", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(args.output_dir / "TOC-Graphic.tif", dpi=args.dpi, facecolor="white",
                bbox_inches="tight", pad_inches=0.02,
                pil_kwargs={"compression": "tiff_lzw"})
    print(f"wrote {args.output_dir}/TOC-Graphic.pdf|tif|png at "
          f"{WIDTH_IN} x {HEIGHT_IN} in, {args.dpi} dpi")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
