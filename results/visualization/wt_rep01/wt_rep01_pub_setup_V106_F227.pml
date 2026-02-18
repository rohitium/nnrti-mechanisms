reinitialize

# Load trajectory (WT rep01)
load /Users/rohitpro/Career/00_Github/nnrti-mechanisms/results/md_runs/wt/rep_01/wt_rep01_analysis_topology.pdb, wt_rep01
load_traj /Users/rohitpro/Career/00_Github/nnrti-mechanisms/results/md_runs/wt/rep_01/wt_rep01_analysis.dcd, wt_rep01, state=1

# WT-equivalent residue/ligand selections for the V106I+F227C context
select v106_res, wt_rep01 and chain A and resn VAL and resi 103
select f227_res, wt_rep01 and chain A and resn PHE and resi 224
select dor, wt_rep01 and resn 2KW

# Sidechain heavy-atom selections for metric-consistent min distances
select v106_sc, v106_res and not name N+CA+C+O+OXT and not elem H
select f227_sc, f227_res and not name N+CA+C+O+OXT and not elem H
select dor_heavy, dor and not elem H

# Label anchors
select v106_anchor, v106_res and name CB
select f227_anchor, f227_res and name CZ
select dor_anchor, first (dor_heavy)

# Display
hide everything
show cartoon, wt_rep01 and polymer.protein
show sticks, dor or v106_res or f227_res
set stick_radius, 0.22, dor or v106_res or f227_res
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_transparency, 0.15

# WT color palette (distinct from mutant)
color gray70, wt_rep01 and polymer.protein
color teal, v106_res
color purple, f227_res
color tv_yellow, dor

# Distance labels are built below as per-frame minimum heavy-atom distances
set dynamic_measures, on
set label_digits, 2

# Persistent labels
set label_font_id, 7
set label_size, 26
set label_color, black
set label_outline_color, white
label v106_anchor, "V106"
label f227_anchor, "F227"
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

# Camera framing around interaction site
orient dor or v106_res or f227_res
zoom dor or v106_res or f227_res, 12
scene pub_site, store

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
    for obj in ["v106_to_dor_min", "f227_to_dor_min"]:
        cmd.delete(obj)

    n_states = cmd.count_states("wt_rep01")
    v106_a_by_state = {}
    v106_b_by_state = {}
    f227_a_by_state = {}
    f227_b_by_state = {}

    for st in range(1, n_states + 1):
        p1 = _min_pair_coords("v106_sc", "dor_heavy", st)
        p2 = _min_pair_coords("f227_sc", "dor_heavy", st)
        if p1 is None or p2 is None:
            continue
        v106_a, v106_b, _ = p1
        f227_a, f227_b, _ = p2
        v106_a_by_state[st] = v106_a
        v106_b_by_state[st] = v106_b
        f227_a_by_state[st] = f227_a
        f227_b_by_state[st] = f227_b

    _build_multistate_point("v106_min_a_pt", v106_a_by_state, n_states)
    _build_multistate_point("v106_min_b_pt", v106_b_by_state, n_states)
    _build_multistate_point("f227_min_a_pt", f227_a_by_state, n_states)
    _build_multistate_point("f227_min_b_pt", f227_b_by_state, n_states)

    cmd.distance("v106_to_dor_min", "v106_min_a_pt", "v106_min_b_pt", state=0)
    cmd.distance("f227_to_dor_min", "f227_min_a_pt", "f227_min_b_pt", state=0)
    cmd.hide("everything", "v106_min_a_pt or v106_min_b_pt or f227_min_a_pt or f227_min_b_pt")
    cmd.show("dashes", "v106_to_dor_min or f227_to_dor_min")
    cmd.show("labels", "v106_to_dor_min or f227_to_dor_min")
    cmd.color("teal", "v106_to_dor_min")
    cmd.color("purple", "f227_to_dor_min")
    cmd.set("dash_width", 3.0, "v106_to_dor_min or f227_to_dor_min")
    cmd.set("dash_gap", 0.25, "v106_to_dor_min or f227_to_dor_min")
    cmd.set("label_size", 18, "v106_to_dor_min or f227_to_dor_min")
    cmd.set("label_color", "black", "v106_to_dor_min or f227_to_dor_min")

def report_frame_distances(state=1):
    state = int(state)
    cmd.frame(state)
    d1 = cmd.get_distance("v106_min_a_pt", "v106_min_b_pt", state=state)
    d2 = cmd.get_distance("f227_min_a_pt", "f227_min_b_pt", state=state)
    if d1 is None or d2 is None:
        print("Missing measurements at state", state)
        return
    print(f"state={state} V106_to_DOR={d1:.3f}A F227_to_DOR={d2:.3f}A")

def goto_matched_from_mutant(mut_state=280):
    # Mapping from mutant V106I+F227C rep01 -> WT rep01 nearest-time states
    mapping = {260: 361, 280: 389, 310: 431}
    mut_state = int(mut_state)
    if mut_state not in mapping:
        print("Known mutant states:", sorted(mapping.keys()))
        return
    wt_state = mapping[mut_state]
    cmd.frame(wt_state)
    print(f"mutant_state={mut_state} -> wt_state={wt_state}")

cmd.extend("build_min_distance_objects", build_min_distance_objects)
cmd.extend("report_frame_distances", report_frame_distances)
cmd.extend("goto_matched_from_mutant", goto_matched_from_mutant)
build_min_distance_objects()
python end

scene pub_site, store
