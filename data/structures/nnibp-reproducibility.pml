# Is the 4NCG doravirine pose a one-off, or a reproducible property of the pocket?
#
# Three renders, all from ONE superposition: the 16 NNIBP-lining Calpha atoms of
# 4NCG (15 in p66, Glu138 in p51). Every structure was written already in that
# frame by nnrti.cli.compare_nnrti_crystal_poses, so the pictures and the
# numbers in Supplementary Tables 8 and 9 come from the same fit.
#
#   A  the two independent doravirine structures   -> the ligand agrees
#   B  Tyr181 across five NNRTI complexes          -> the pocket CAN move
#   C  Tyr188 across the same five                 -> and here it does not
#
# B is the yardstick that makes A meaningful: the same residue, under the same
# superposition, adopts a different rotamer with a different drug.
#
# Run from the repository root:
#   /Applications/PyMOL.app/Contents/MacOS/PyMOL -cq data/structures/nnibp-reproducibility.pml

reinitialize
bg_color white
set ray_opaque_background, 0
set opaque_background, 0
set ray_shadows, 0
set antialias, 2
set two_sided_lighting, on
set ray_interior_color, grey90
set transparency_mode, 2
set cartoon_discrete_colors, 1
set cartoon_transparency, 0.9
set cartoon_side_chain_helper, 1
set valence, 0
set stick_radius, 0.20
set stick_h_scale, 1.0

python
from pymol import cmd

D = "results/analysis/nnrti_pose_comparison/structures"
OUT = "results/analysis/nnrti_pose_comparison/plots"
NNIBP = "100+101+103+106+107+108+179+181+188+189+190+227+229+234+318"

STRUCTS = [
    ("dor_xray", "4NCG_on_4NCG_nnibp_frame.pdb", "2KW", "grey50"),
    ("dor_em",   "7Z2G_on_4NCG_nnibp_frame.pdb", "2KW", "marine"),
    ("rpv",      "2ZD1_on_4NCG_nnibp_frame.pdb", "T27", "firebrick"),
    ("efv",      "1FK9_on_4NCG_nnibp_frame.pdb", "EFZ", "orange"),
    ("nvp",      "1VRT_on_4NCG_nnibp_frame.pdb", "NVP", "purple"),
]

for obj, f, lig, col in STRUCTS:
    cmd.load(f"{D}/{f}", obj)
    cmd.remove(f"{obj} and solvent")
    cmd.remove(f"{obj} and hydro")
    cmd.remove(f"{obj} and polymer.nucleic")       # 7Z2G carries a DNA aptamer
    cmd.hide("everything", obj)
    cmd.select(f"{obj}_lig", f"{obj} and resn {lig}")
    cmd.color(col, f"{obj}_lig and elem C")
    cmd.util.cnc(f"{obj}_lig")


def shot(name, w=1800, h=1500):
    cmd.png(f"{OUT}/{name}", width=w, height=h, dpi=300, ray=1)


# ---------------------------------------------------------------- A --------
# the two doravirine structures, in the pocket that holds them
cmd.hide("everything")
cmd.show("cartoon", "dor_xray and chain A")
cmd.color("grey90", "dor_xray and chain A")
cmd.show("sticks", f"dor_xray and chain A and resi {NNIBP} and sidechain")
cmd.set("stick_radius", 0.10, "dor_xray")
cmd.color("grey75", f"dor_xray and chain A and resi {NNIBP} and elem C")
cmd.set("stick_transparency", 0.45, f"dor_xray and chain A and resi {NNIBP}")
for o in ("dor_xray", "dor_em"):
    cmd.show("sticks", f"{o}_lig")
    cmd.set("stick_radius", 0.22, f"{o}_lig")
cmd.orient("dor_xray_lig")
cmd.turn("y", 15)
cmd.zoom("dor_xray_lig", 4.0)
shot("A_two_DOR_structures.png")


# --------------------------------------------------------------- B, C ------
# one residue at a time, nothing else in the scene
def residue_panel(resi, name, turn_y=0):
    cmd.hide("everything")
    for obj, _f, _l, col in STRUCTS:
        sel = f"{obj} and chain A and resi {resi} and sidechain"
        cmd.show("sticks", sel)
        cmd.set("stick_radius", 0.24, sel)
        cmd.color(col, f"{sel} and elem C")
        cmd.util.cnc(sel)
    # the doravirine molecule, faint, for orientation only
    cmd.show("sticks", "dor_xray_lig")
    cmd.set("stick_radius", 0.11, "dor_xray_lig")
    cmd.color("grey80", "dor_xray_lig and elem C")
    cmd.set("stick_transparency", 0.55, "dor_xray_lig")
    cmd.orient(f"dor_xray and chain A and resi {resi}")
    if turn_y:
        cmd.turn("y", turn_y)
    cmd.zoom(f"dor_xray and chain A and resi {resi}", 3.2)
    shot(name)


residue_panel(181, "B_Tyr181_rotamer.png", turn_y=25)
residue_panel(188, "C_Tyr188_rotamer.png", turn_y=25)
print("wrote A, B, C")
python end
