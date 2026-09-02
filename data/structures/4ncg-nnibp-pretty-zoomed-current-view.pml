# Apply publication styling to the CURRENT 4ncg view without resetting camera.
# This preserves your hand-tuned viewpoint and only zooms in.

python
from pymol import cmd

if cmd.count_atoms("4ncg") == 0:
    cmd.load("data/structures/4NCG.cif", "4ncg")
python end

# Core selections
select dor, 4ncg and resn 2KW
select dor_heavy, dor and not elem H
select key_rt, 4ncg and chain A and resi 100+103+106+188+227+229+234+236
select key_rt_ca, key_rt and name CA

# Clean previous helper/measurement objects
python
from pymol import cmd
for name in cmd.get_names("objects"):
    if name.startswith("dist_") or name.endswith("_min_a") or name.endswith("_min_b") or name == "dor_label_pt":
        cmd.delete(name)
python end

hide everything, 4ncg
show sticks, dor
show lines, key_rt
set stick_radius, 0.24, dor
set line_width, 3.0, key_rt

# Element coloring for BOTH residues and DOR
color green, (dor or key_rt) and elem C
color blue, (dor or key_rt) and elem N
color red, (dor or key_rt) and elem O
color yellow, (dor or key_rt) and elem S
color cyan, (dor or key_rt) and elem F+CL
color orange, (dor or key_rt) and elem P

# Labels
set label_font_id, 7
set label_size, 22
set label_color, black
set label_outline_color, white
set label_digits, 1
set float_labels, on
label key_rt_ca, "%s%s" % (resn,resi)

# Publication render settings
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

# DOR label anchor
com_dor = cmd.centerofmass("dor_heavy")
com_key = cmd.centerofmass("key_rt_ca")
vec = cpv.sub(com_dor, com_key)
if cpv.length(vec) < 1.0e-6:
    vec = [0.0, 0.0, 1.0]
vec = cpv.normalize(vec)
pos = cpv.add(com_dor, cpv.scale(vec, 1.8))
cmd.pseudoatom("dor_label_pt", pos=list(pos))
cmd.hide("everything", "dor_label_pt")
cmd.label("dor_label_pt", '"DOR"')

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
    base = f"{resn}{resi}"
    a_pt = f"{base}_min_a"
    b_pt = f"{base}_min_b"
    d_obj = f"dist_{base}_DOR"

    cmd.pseudoatom(a_pt, pos=list(a))
    cmd.pseudoatom(b_pt, pos=list(b))
    cmd.distance(d_obj, a_pt, b_pt)
    cmd.hide("everything", f"{a_pt} or {b_pt}")

cmd.show("dashes", "dist_*")
cmd.show("labels", "dist_*")
cmd.color("yellow", "dist_*")
cmd.set("dash_width", 3.0, "dist_*")
cmd.set("dash_gap", 0.20, "dist_*")
cmd.set("label_size", 18, "dist_*")
cmd.set("label_color", "black", "dist_*")
python end

# Keep current camera orientation, but zoom in tighter
zoom dor or key_rt, 8

# High-res transparent export
ray 3600, 2700
png data/structures/4ncg-nnibp-pretty-zoomed-current-view.png, dpi=600
