reinitialize

# Load trajectory (V106I+F227C rep01)
load /Users/rohitpro/Career/00_Github/nnrti-mechanisms/results/md_runs/V106I_F227C/rep_01/V106I_F227C_rep01_analysis_topology.pdb, V106I_F227C_rep01
load_traj /Users/rohitpro/Career/00_Github/nnrti-mechanisms/results/md_runs/V106I_F227C/rep_01/V106I_F227C_rep01_analysis.dcd, V106I_F227C_rep01, state=1

# Key residue/ligand selections
select v106i_res, V106I_F227C_rep01 and chain A and resn ILE and resi 103
select f227c_res, V106I_F227C_rep01 and chain A and resn CYS and resi 224
select dor, V106I_F227C_rep01 and chain C and resn 2KW

# Sidechain selections (heavy atoms only) for contact context
select v106i_sc, v106i_res and not name N+CA+C+O+OXT and not elem H
select f227c_sc, f227c_res and not name N+CA+C+O+OXT and not elem H
select dor_heavy, dor and not elem H

# Label anchors
select v106i_anchor, v106i_res and name CB
select f227c_anchor, f227c_res and name SG
select dor_anchor, first (dor_heavy)

# Display
hide everything
show cartoon, V106I_F227C_rep01 and polymer.protein
show sticks, dor or v106i_res or f227c_res
set stick_radius, 0.22, dor or v106i_res or f227c_res
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_transparency, 0.15

# Colors
color gray80, V106I_F227C_rep01 and polymer.protein
color tv_orange, v106i_res
color marine, f227c_res
color yelloworange, dor

# Per-frame minimum sidechain-heavy-atom to DOR-heavy-atom distances are
# generated in the Python block below (objects: v106i_to_dor_min, f227c_to_dor_min).
set dynamic_measures, on
set label_digits, 2

# Persistent labels
set label_font_id, 7
set label_size, 26
set label_color, black
set label_outline_color, white
label v106i_anchor, "V106I"
label f227c_anchor, "F227C"
label dor_anchor, "DOR"

# Publication-style rendering
set orthoscopic, on
set depth_cue, 0
set antialias, 2
set ray_opaque_background, off
set ray_shadow, 0
set ray_trace_mode, 1
set ray_trace_gain, 0.08
set ambient, 0.2
set direct, 0.75
set spec_reflect, 0.25
set spec_power, 80
bg_color white

# Camera framing around the interaction site
orient dor or v106i_res or f227c_res
zoom dor or v106i_res or f227c_res, 12
scene pub_site, store

# Optional helper commands (no files written unless you call png/save manually)
python
from pymol import cmd
from chempy import cpv

def _min_pair_coords(sel_a, sel_b, state):
    a_atoms = cmd.get_model(sel_a, state).atom
    b_atoms = cmd.get_model(sel_b, state).atom
    if not a_atoms or not b_atoms:
        return None
    best_a = None
    best_b = None
    best_d = 1.0e9
    for a in a_atoms:
        ac = a.coord
        for b in b_atoms:
            d = cpv.distance(ac, b.coord)
            if d < best_d:
                best_d = d
                best_a = ac
                best_b = b.coord
    return best_a, best_b, best_d

def _build_multistate_point(obj_name, coords_by_state, n_states):
    cmd.delete(obj_name)
    first = coords_by_state.get(1)
    if first is None:
        first = next(iter(coords_by_state.values()))
    cmd.pseudoatom(obj_name, pos=list(first), state=1)
    for st in range(2, n_states + 1):
        cmd.create(obj_name, obj_name, 1, st)
        coord = coords_by_state.get(st, first)
        cmd.alter_state(
            st,
            obj_name,
            "x=coord[0]; y=coord[1]; z=coord[2]",
            space={"coord": coord},
        )

def build_min_distance_objects():
    for obj in ["v106i_to_dor_min", "f227c_to_dor_min"]:
        cmd.delete(obj)

    n_states = cmd.count_states("V106I_F227C_rep01")
    v106i_a_by_state = {}
    v106i_b_by_state = {}
    f227c_a_by_state = {}
    f227c_b_by_state = {}

    for st in range(1, n_states + 1):
        p1 = _min_pair_coords("v106i_sc", "dor_heavy", st)
        p2 = _min_pair_coords("f227c_sc", "dor_heavy", st)
        if p1 is None or p2 is None:
            continue
        v106i_a, v106i_b, _ = p1
        f227c_a, f227c_b, _ = p2
        v106i_a_by_state[st] = v106i_a
        v106i_b_by_state[st] = v106i_b
        f227c_a_by_state[st] = f227c_a
        f227c_b_by_state[st] = f227c_b

    _build_multistate_point("v106i_min_a_pt", v106i_a_by_state, n_states)
    _build_multistate_point("v106i_min_b_pt", v106i_b_by_state, n_states)
    _build_multistate_point("f227c_min_a_pt", f227c_a_by_state, n_states)
    _build_multistate_point("f227c_min_b_pt", f227c_b_by_state, n_states)

    cmd.distance("v106i_to_dor_min", "v106i_min_a_pt", "v106i_min_b_pt", state=0)
    cmd.distance("f227c_to_dor_min", "f227c_min_a_pt", "f227c_min_b_pt", state=0)
    cmd.hide("everything", "v106i_min_a_pt or v106i_min_b_pt or f227c_min_a_pt or f227c_min_b_pt")
    cmd.show("dashes", "v106i_to_dor_min or f227c_to_dor_min")
    cmd.show("labels", "v106i_to_dor_min or f227c_to_dor_min")
    cmd.color("tv_orange", "v106i_to_dor_min")
    cmd.color("marine", "f227c_to_dor_min")
    cmd.set("dash_width", 3.0, "v106i_to_dor_min or f227c_to_dor_min")
    cmd.set("dash_gap", 0.25, "v106i_to_dor_min or f227c_to_dor_min")
    cmd.set("label_size", 18, "v106i_to_dor_min or f227c_to_dor_min")
    cmd.set("label_color", "black", "v106i_to_dor_min or f227c_to_dor_min")

def report_frame_distances(state=1):
    state = int(state)
    cmd.frame(state)
    d1 = cmd.get_distance("v106i_min_a_pt", "v106i_min_b_pt", state=state)
    d2 = cmd.get_distance("f227c_min_a_pt", "f227c_min_b_pt", state=state)
    if d1 is None or d2 is None:
        print("Missing measurements at state", state)
        return
    print(f"state={state} V106I_to_DOR={d1:.3f}A F227C_to_DOR={d2:.3f}A")

cmd.extend("build_min_distance_objects", build_min_distance_objects)
cmd.extend("report_frame_distances", report_frame_distances)
build_min_distance_objects()
python end

# Re-store scene after creating min-distance objects
scene pub_site, store
