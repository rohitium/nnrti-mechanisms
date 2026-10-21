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


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pml", type=Path, default=Path("ops/figures/toc_graphic.pml"))
    p.add_argument("--scene", type=Path, default=Path("/tmp/toc_scene.png"))
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

    ax = fig.add_axes([0.0, 0.0, 0.545, 1.0])
    ax.imshow(img, interpolation="lanczos")
    ax.axis("off")

    # Text lives in its own axes in data coordinates. Figure-fraction placement
    # stretches a circle into an ellipse and makes it easy to run a title off
    # the canvas, both of which happened in the first attempt.
    tb = fig.add_axes([0.555, 0.0, 0.445, 1.0])
    tb.set_xlim(0, 1)
    tb.set_ylim(0, 1)
    tb.axis("off")

    tb.text(0.0, 0.93, "Doravirine resistance in\nHIV-1 reverse transcriptase",
            fontsize=7.6, fontweight="bold", color="#141414", va="top", linespacing=1.3)
    tb.text(0.0, 0.645,
            "18 clinical genotypes scored by\nalchemical free energy perturbation\n"
            "and MM/GBSA, 3 \u00d7 100 ns each",
            fontsize=5.8, color="#3d3d3d", va="top", linespacing=1.65)
    tb.scatter([0.028], [0.215], s=48, color="#f08c1b", edgecolor="#6b3300",
               linewidths=0.5, zorder=5)
    tb.text(0.095, 0.215, "substituted positions", fontsize=5.8, color="#3d3d3d",
            va="center")
    tb.text(0.0, 0.075, "doravirine bound in the NNIBP", fontsize=5.6,
            color="#6a6a6a", va="center")

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
