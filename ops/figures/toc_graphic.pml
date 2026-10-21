# TOC graphic for ci-2026-030865: doravirine bound in the NNIBP of HIV-1 RT.
#
# Renders the molecular scene only. Text is composited afterwards at known
# physical size by plot_toc_graphic.py, so labels stay legible at the 3.25 in
# printed width instead of being baked in at whatever scale the ray happened
# to use.
#
# Run via: PYTHONPATH=src python -m nnrti.cli.plot_toc_graphic

load results/md_runs/wt/rep_01/wt_minimized_rep01.pdb, sys
hide everything
remove solvent
remove not polymer and not resn 2KW
# This is a minimised MD structure, so every hydrogen is explicit. Left in,
# the side chains render as a thicket of lines that buries the ligand.
remove hydrogens

select p66, polymer and chain A
select dor, resn 2KW
select lining, byres (p66 within 5.0 of dor)

# Substituted positions, topology numbering (canonical - 3). Split into their
# own object so colouring the spheres cannot bleed into the cartoon ribbon.
create sites, p66 and name CA and resi 95+97+100+103+178+185+187+222+224+227+231+315

bg_color white
set ray_opaque_background, 1
set antialias, 2
set ray_trace_mode, 1
set ray_trace_gain, 0.10
set ray_trace_color, grey30
set ray_shadows, 0
set depth_cue, 0.35
set fog_start, 0.60
set specular, 0.30
set surface_quality, 1

# Context ribbon, deliberately pale: the cleft should read, not the fold.
show cartoon, p66
color grey85, p66
set cartoon_transparency, 0.52, p66

# Only the residues the paper actually argues about, as thin sticks. Showing
# every pocket-lining side chain buried the ligand in line noise.
select keyres, p66 and resi 100+103+178+185+187 and sidechain
show sticks, keyres
color palecyan, keyres and elem C
set stick_radius, 0.085, keyres

# Doravirine: the subject of the figure, so it gets the strongest contrast.
show sticks, dor
set stick_radius, 0.21, dor
color yelloworange, dor and elem C
color tv_blue, dor and elem N
color firebrick, dor and elem O
color limon, dor and elem F
color green, dor and elem Cl

show spheres, sites
set sphere_scale, 0.40, sites
color orange, sites
set sphere_transparency, 0.0, sites

orient dor
zoom dor, 10
turn x, -6
turn y, 20
set field_of_view, 20

ray 1180, 1050
png /tmp/toc_scene.png, dpi=600
