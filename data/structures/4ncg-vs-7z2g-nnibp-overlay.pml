# Overlay of the two experimental RT-DOR structures on the NNIBP.
#
#   4NCG  2.58 A X-ray,   RT + DOR
#   7Z2G  3.65 A cryo-EM, RT + DOR + 38-mer DNA aptamer
#
# The fit is on the 16 NNIBP-lining Ca atoms only (15 p66 + Glu138 of p51), so
# what the picture shows is whether DOR occupies the same place relative to the
# residues that line it -- not whether the two enzymes superpose globally, which
# they do not, because only one of them has nucleic acid bound.
#
# Numbers behind the picture: results/analysis/crystal_pose_comparison/
#
# Run from the repository root:
#   pymol -cq data/structures/4ncg-vs-7z2g-nnibp-overlay.pml

reinitialize

load data/structures/apo_reference/4NCG.pdb, xray
# already placed in the NNIBP frame by nnrti.cli.compare_dor_crystal_poses, so
# the picture and the numbers come from one and the same superposition
load results/analysis/crystal_pose_comparison/structures/7Z2G_on_4NCG_nnibp_frame.pdb, cryoem

remove solvent
remove hydrogens

# NNIBP lining, canonical HIV-1 RT numbering (Cilento, Kirby & Sarafianos)
select nnibp_xray,   xray   and chain A and resi 100+101+103+106+107+108+179+181+188+189+190+227+229+234+318
select nnibp_cryoem, cryoem and chain A and resi 100+101+103+106+107+108+179+181+188+189+190+227+229+234+318

select dor_xray,   xray   and resn 2KW
select dor_cryoem, cryoem and resn 2KW

hide everything
show sticks, dor_xray or dor_cryoem
show lines, (nnibp_xray or nnibp_cryoem) and sidechain
show cartoon, xray and chain A and resi 95-115+175-195+225-240+315-320

# the ligands carry the comparison -- the pocket is context, drawn thin
set cartoon_transparency, 0.9
set stick_radius, 0.22, dor_xray or dor_cryoem
set line_width, 1.6
set valence, 0

color grey90, xray
color grey90, cryoem
color grey70, (nnibp_xray or nnibp_cryoem) and elem C
color marine,    dor_xray and elem C
color firebrick, dor_cryoem and elem C
util.cnc("(dor_xray or dor_cryoem) and not elem C")

set orthoscopic, on
set ray_opaque_background, 1
set ray_shadows, 0
set antialias, 2
set depth_cue, 0
bg_color white

orient dor_xray
zoom dor_xray, 3.0
turn x, -10

png results/analysis/crystal_pose_comparison/plots/nnibp_overlay.png, width=2400, height=1800, dpi=300, ray=1

# a second view down the long axis of the ligand, where a tilt would show up
turn y, 90
png results/analysis/crystal_pose_comparison/plots/nnibp_overlay_side.png, width=2400, height=1800, dpi=300, ray=1
