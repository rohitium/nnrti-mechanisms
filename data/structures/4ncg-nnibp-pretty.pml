# Run from the repository root:  pymol -cq data/structures/4ncg-nnibp-pretty.pml
reinitialize

# Load crystal structure (4NCG, Doravirine = resn 2KW)
load data/structures/4NCG.cif, 4ncg

# Core selections
select dor, 4ncg and resn 2KW
select dor_heavy, dor and not elem H
select key_rt, 4ncg and chain A and resi 100+103+106+188+227+229+234+236

# Label anchors
select key_rt_ca, key_rt and name CA

# Display cleanup
hide everything
show sticks, dor
show lines, key_rt

# Styling
set stick_radius, 0.24, dor
set line_width, 3.0, key_rt
color tv_orange, dor
color forest, key_rt

# Label styling for readability
set label_font_id, 7
set label_size, 24
set label_color, black
set label_outline_color, white
set label_bg_color, white
set label_bg_transparency, 0.55
set label_digits, 2
set float_labels, on

# Static residue labels
label key_rt_ca, "%s%s" % (resn,resi)

# Render settings (publication-style, transparent background)
set orthoscopic, on
set depth_cue, 0
set antialias, 2
set ray_opaque_background, off
set ray_shadow, 0
set ray_trace_mode, 1
set ray_trace_gain, 0.08
set ambient, 0.22
set direct, 0.75
set spec_reflect, 0.22
set spec_power, 80
bg_color white

# Place a clean DOR label on a dedicated pseudoatom to avoid overlap
python
from pymol import cmd
from chempy import cpv

def _min_pair_coords(sel_a, sel_b, state=1):
    a_atoms = cmd.get_model(sel_a, state).atom
    b_atoms = cmd.get_model(sel_b, state).atom
    if not a_atoms or not b_atoms:
        return None
    best = None
    best_d = 1.0e9
    for a in a_atoms:
        ac = a.coord
        for b in b_atoms:
            d = cpv.distance(ac, b.coord)
            if d < best_d:
                best_d = d
                best = (ac, b.coord, d)
    return best

# DOR label anchor (offset from DOR center toward camera-facing side)
cmd.delete("dor_label_pt")
com_dor = cmd.centerofmass("dor_heavy")
com_key = cmd.centerofmass("key_rt_ca")
vec = cpv.sub(com_dor, com_key)
if cpv.length(vec) < 1.0e-6:
    vec = [0.0, 0.0, 1.0]
vec = cpv.normalize(vec)
pos = cpv.add(com_dor, cpv.scale(vec, 2.2))
cmd.pseudoatom("dor_label_pt", pos=list(pos))
cmd.hide("everything", "dor_label_pt")
cmd.label("dor_label_pt", '"DOR"')

# Build minimum all-heavy-atom (backbone + sidechain) distance object per key residue
for name in cmd.get_names("objects"):
    if name.startswith("dist_") or name.endswith("_min_a") or name.endswith("_min_b"):
        cmd.delete(name)

key_residues = [
    ("LEU", "100"),
    ("LYS", "103"),
    ("VAL", "106"),
    ("TYR", "188"),
    ("PHE", "227"),
    ("TRP", "229"),
    ("LEU", "234"),
    ("PRO", "236"),
]

for resn, resi in key_residues:
    res_sel = f"4ncg and chain A and resi {resi} and resn {resn}"
    residue_heavy = f"({res_sel}) and not elem H"

    pair = _min_pair_coords(residue_heavy, "dor_heavy", state=1)
    if pair is None:
        continue

    a, b, _ = pair
    obj_base = f"{resn}{resi}"
    a_pt = f"{obj_base}_min_a"
    b_pt = f"{obj_base}_min_b"
    d_obj = f"dist_{obj_base}_DOR"

    cmd.pseudoatom(a_pt, pos=list(a))
    cmd.pseudoatom(b_pt, pos=list(b))
    cmd.distance(d_obj, a_pt, b_pt)
    cmd.hide("everything", f"{a_pt} or {b_pt}")

# Distance object styling
cmd.show("dashes", "dist_*")
cmd.show("labels", "dist_*")
cmd.color("yellow", "dist_*")
cmd.set("dash_width", 3.0, "dist_*")
cmd.set("dash_gap", 0.20, "dist_*")
cmd.set("label_size", 18, "dist_*")
cmd.set("label_color", "black", "dist_*")
python end

# Camera framing
orient dor or key_rt
zoom dor or key_rt, 12

# Save a publication PNG with transparent background
ray 3600, 2700
png data/structures/4ncg-nnibp-pretty.png, dpi=600
