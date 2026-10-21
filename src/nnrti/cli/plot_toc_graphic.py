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
#: ACS gives 3.25 x 1.75 in as the maximum TOC size. The height here is
#: derived from the panel band rather than fixed, so it comes out under it.
WIDTH_IN = 3.25


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

    cycle = mpimg.imread(str(trim_white(args.cycle)))

    # Both panels share one vertical band, so their tops and bottoms line up.
    # Centring each independently leaves mismatched margins, because the
    # molecular render is square and the cycle is about 1.5:1 wide. The canvas
    # height is then derived from the content rather than fixed in advance: a
    # fixed height taller than the band would simply reinstate the whitespace
    # this is removing. ACS gives a maximum TOC size, not a required one.
    a_mol = img.shape[1] / img.shape[0]
    a_cyc = cycle.shape[1] / cycle.shape[0]

    margin, gap = 0.022, 0.018                       # figure fractions of width
    usable_w_in = WIDTH_IN * (1.0 - 2 * margin - gap)
    band_h_in = usable_w_in / (a_mol + a_cyc)        # common panel height
    height_in = band_h_in + 2 * (margin * WIDTH_IN)  # equal margins top/bottom

    fig = plt.figure(figsize=(WIDTH_IN, height_in))
    fig.patch.set_facecolor("white")

    band_h = band_h_in / height_in
    y0 = (1.0 - band_h) / 2.0
    w_mol = (band_h_in * a_mol) / WIDTH_IN
    w_cyc = (band_h_in * a_cyc) / WIDTH_IN

    ax = fig.add_axes([margin, y0, w_mol, band_h])
    ax.imshow(img, interpolation="lanczos")
    ax.axis("off")

    ax2 = fig.add_axes([margin + w_mol + gap, y0, w_cyc, band_h])
    ax2.imshow(cycle, interpolation="lanczos")
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
          f"{WIDTH_IN:.2f} x {height_in:.2f} in, {args.dpi} dpi "
          f"(panel band {band_h_in:.2f} in, both panels aligned)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
