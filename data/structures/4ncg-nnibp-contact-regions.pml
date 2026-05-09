reinitialize

# Crystal NNIBP view keyed to the WT all-contacted-residue plot.
# Doravirine is resn 2KW in 4NCG. Contacted residues are grouped by the
# same seven plot regions used in plot_triplet_contact_story.py.
load /Users/rohitpro/Career/00_Github/nnrti-mechanisms/data/structures/4NCG.cif, dor_rt

hide everything
remove solvent

# Core selections
select dor, dor_rt and resn 2KW
select dor_heavy, dor and not elem H
select rt_chain_a, dor_rt and chain A and polymer.protein
select pocket_context, rt_chain_a and byres (rt_chain_a within 14 of dor_heavy)

select region_beta6_strand, rt_chain_a and resi 95+97
select region_pocket_entrance, rt_chain_a and resi 100+101+103+179+181
select region_loop_103_108, rt_chain_a and resi 102+104+105+106+107+108
select region_hydrophobic_tunnel, rt_chain_a and resi 188+227+229+234
select region_beta9_beta10_hairpin, rt_chain_a and resi 180+189+190
select region_primer_grip, rt_chain_a and resi 223+225+228+235+236+237
select region_distal_wall, rt_chain_a and resi 318
select contact_regions, region_beta6_strand or region_pocket_entrance or region_loop_103_108 or region_hydrophobic_tunnel or region_beta9_beta10_hairpin or region_primer_grip or region_distal_wall

# Muted pocket context so the region surfaces define the visual story.
show cartoon, pocket_context
color grey80, pocket_context
set cartoon_transparency, 0.95, pocket_context

# Region sidechains and translucent contact-region surfaces.
show sticks, contact_regions and not name N+C+O
set stick_radius, 0.10, contact_regions
show surface, contact_regions
set transparency, 0.50, contact_regions
set surface_quality, 1
set solvent_radius, 1.4

set_color c_beta6, [0.30, 0.47, 0.66]
set_color c_entrance, [0.96, 0.49, 0.11]
set_color c_loop, [0.93, 0.78, 0.21]
set_color c_tunnel, [0.33, 0.62, 0.31]
set_color c_hairpin, [0.89, 0.34, 0.34]
set_color c_primer, [0.69, 0.48, 0.63]
set_color c_distal, [0.45, 0.72, 0.70]

color c_beta6, region_beta6_strand
color c_entrance, region_pocket_entrance
color c_loop, region_loop_103_108
color c_tunnel, region_hydrophobic_tunnel
color c_hairpin, region_beta9_beta10_hairpin
color c_primer, region_primer_grip
color c_distal, region_distal_wall

# DOR is the anchor: slightly larger sticks plus atom-type colors.
show sticks, dor
hide spheres, dor
set stick_radius, 0.34, dor
color yellow, dor and elem C
color blue, dor and elem N
color red, dor and elem O
color cyan, dor and elem F
color green, dor and elem Cl

# Crisp outlines and lighting.
set orthoscopic, on
set antialias, 2
set depth_cue, 0
set ambient, 0.28
set direct, 0.72
set spec_reflect, 0.22
set spec_power, 80
set ray_opaque_background, off
set ray_shadow, 0
set ray_trace_mode, 0
set ray_trace_gain, 0.04
set two_sided_lighting, on
bg_color white

# Large dark region labels, offset away from DOR so they annotate the colored surfaces.
set label_font_id, 7
set label_size, 15
set label_color, black
set label_outline_color, black
set float_labels, on

python
from pymol import cmd
from chempy import cpv

for name in list(cmd.get_names("objects")):
    if name.startswith("label_"):
        cmd.delete(name)

def label_region(name, selection, text, scale=4.2, manual=(0.0, 0.0, 0.0)):
    if cmd.count_atoms(selection) == 0:
        return
    region_com = cmd.centerofmass(selection)
    dor_com = cmd.centerofmass("dor_heavy")
    direction = cpv.sub(region_com, dor_com)
    if cpv.length(direction) < 1.0e-6:
        direction = [0.0, 0.0, 1.0]
    direction = cpv.normalize(direction)
    pos = cpv.add(region_com, cpv.add(cpv.scale(direction, scale), list(manual)))
    obj = "label_" + name
    cmd.pseudoatom(obj, pos=list(pos))
    cmd.hide("everything", obj)
    cmd.label(obj, repr(text))

label_region("beta6", "region_beta6_strand", "β6 strand", 4.0, (-1.4, 0.2, 0.8))
label_region("entrance", "region_pocket_entrance", "Pocket\nentrance", 4.4, (-1.0, -0.9, 0.5))
label_region("loop103108", "region_loop_103_108", "103-108 loop", 4.3, (0.4, -1.2, 1.0))
label_region("tunnel", "region_hydrophobic_tunnel", "Hydrophobic\ntunnel", 4.6, (1.0, 0.7, 0.8))
label_region("hairpin", "region_beta9_beta10_hairpin", "β9-β10\nhairpin", 4.5, (-1.8, 1.5, 0.8))
label_region("primer", "region_primer_grip", "β12-β13\nprimer grip", 5.0, (1.2, 1.5, 0.8))
label_region("distal", "region_distal_wall", "Distal\npocket wall", 4.8, (0.6, -0.6, 1.2))
python end

# Camera: frame the ligand lengthwise inside the region surfaces.
orient dor or contact_regions
turn x, 12
turn y, -32
turn z, 48
zoom dor or contact_regions, 14
move x, -8
move y, 1
clip slab, 42

ray 3600, 2200
png /Users/rohitpro/Career/00_Github/nnrti-mechanisms/data/structures/4ncg-nnibp-contact-regions.png, dpi=600
