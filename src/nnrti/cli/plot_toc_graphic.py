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


def panel_a(pdf: Path, dpi: int = 150) -> Path:
    """Crop panel A, the thermodynamic cycle, out of Figure 2.

    The whole four-panel figure cannot go in a TOC graphic: it is 13 in wide in
    the paper and the TOC is 3.3 in, so its 8 pt labels would land near 1 pt and
    fail the legibility requirement. Panel A alone carries the point -- the
    alchemical cycle that makes ddG_bind an equilibrium quantity -- at a size
    that reads.
    """
    from PIL import Image
    out = Path("/tmp/toc_fig2_panelA.png")
    full = Path("/tmp/toc_fig2_full.png")
    subprocess.run(["pdftoppm", "-png", "-r", str(dpi), "-singlefile",
                    str(pdf), str(full.with_suffix(""))], check=True)
    im = Image.open(full)
    w, h = im.size
    im.crop((int(0.118 * w), int(0.020 * h), int(0.503 * w), int(0.492 * h))).save(out)
    return out


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pml", type=Path, default=Path("ops/figures/toc_graphic.pml"))
    p.add_argument("--scene", type=Path, default=Path("/tmp/toc_scene.png"))
    p.add_argument("--figure2-pdf", type=Path,
                   default=Path("paper/submission/figures/Fig2.pdf"),
                   help="Figure 2 of the manuscript; panel A is cropped from it.")
    p.add_argument("--figure2", type=Path, default=None,
                   help="Use this image for the right-hand panel instead of cropping "
                        "panel A from the Figure 2 PDF.")
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

    fig2 = mpimg.imread(str(args.figure2 or panel_a(args.figure2_pdf)))
    h, w = fig2.shape[:2]
    box_w = 0.555
    box_h = box_w * (h / w) * (WIDTH_IN / HEIGHT_IN)
    ax2 = fig.add_axes([0.44, 0.575 - box_h / 2, box_w, box_h])
    ax2.imshow(fig2, interpolation="lanczos")
    ax2.axis("off")

    # The one label on the graphic, bottom right and wrapped so it sits clear
    # of both the structure and the cycle. Its own axes keeps the marker round;
    # a patch placed in figure coordinates gets stretched by the aspect ratio.
    lg = fig.add_axes([0.735, 0.012, 0.255, 0.14])
    lg.set_xlim(0, 1)
    lg.set_ylim(0, 1)
    lg.axis("off")
    lg.scatter([0.08], [0.55], s=28, color="#f08c1b", edgecolor="#6b3300",
               linewidths=0.4, zorder=5)
    lg.text(0.20, 0.55, "substituted\npositions", fontsize=5.2, color="#3d3d3d",
            va="center", ha="left", linespacing=1.35)

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
