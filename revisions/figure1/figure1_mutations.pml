# Figure 1: mutations mapped onto RT-DOR.
# Built by nnrti.cli.build_figure1_pymol -- re-run it to regenerate.
# Everything below is a starting point; adjust and re-save the session.

reinitialize
load /Users/rohitpro/Career/00_Github/nnrti-mechanisms/results/md_runs/wt/rep_01/wt_rep01_analysis_topology.pdb, RT_DOR
hide everything

# ---- the three molecules -------------------------------------------
select p66, RT_DOR and chain A and polymer
select p51, RT_DOR and chain B and polymer
select DOR, RT_DOR and chain E

# ---- the pocket ------------------------------------------------------
select NNIBP, (p66 and resi 97+98+100+103+104+105+176+178+185+186+187+224+226+231+315) or (p51 and resi 135)

# ---- mutated positions ----------------------------------------------
select mutations, p66 and resi 95+97+100+103+178+185+187+222+224+227+231+315
select mut_A98, p66 and resi 95
select mut_L100, p66 and resi 97
select mut_K103, p66 and resi 100
select mut_V106, p66 and resi 103
select mut_Y181, p66 and resi 178
select mut_Y188, p66 and resi 185
select mut_G190, p66 and resi 187
select mut_P225, p66 and resi 222
select mut_F227, p66 and resi 224
select mut_M230, p66 and resi 227
select mut_L234, p66 and resi 231
select mut_Y318, p66 and resi 315

# ---- representations -------------------------------------------------
show cartoon, p66 or p51
color grey80, p66
color grey70, p51
set cartoon_transparency, 0.82, p66
set cartoon_transparency, 0.90, p51
set cartoon_transparency, 0.45, NNIBP

hide everything, DOR and hydro
show sticks, DOR
color yellow,  DOR and elem C
color blue,    DOR and elem N
color red,     DOR and elem O
color aquamarine, DOR and elem F
color green,   DOR and elem Cl
show spheres, DOR
set sphere_scale, 0.20, DOR

show sticks, mutations and not hydro and not (name C+O) and not (name N and not resn PRO)
show spheres, mutations and name CA
set sphere_scale, 0.26, mutations and name CA

show cartoon, NNIBP
color palecyan, NNIBP

# ---- colour the mutated sites --------------------------------------
color orange, mutations
color orange, mutations and elem C
color blue,   mutations and elem N
color red,    mutations and elem O
color yellow, mutations and elem S

# b-factor at each mutated position is log10(fold reduction in DOR
# susceptibility), left in place so the sites can be coloured by
# phenotype instead of flat orange. Uncomment the spectrum line.
alter RT_DOR, b=0.0
alter mut_A98, b=1.968
alter mut_L100, b=0.833
alter mut_K103, b=1.556
alter mut_V106, b=2.207
alter mut_Y181, b=0.146
alter mut_Y188, b=2.173
alter mut_G190, b=1.255
alter mut_P225, b=2.185
alter mut_F227, b=2.025
alter mut_M230, b=1.556
alter mut_L234, b=2.207
alter mut_Y318, b=1.041
rebuild
# spectrum b, yellow_orange_red, mutations, minimum=0.00, maximum=2.21
# ---- labels ----------------------------------------------------------
# Canonical numbering, which is what the manuscript uses.
label mut_A98 and name CA, "A98"
label mut_L100 and name CA, "L100"
label mut_K103 and name CA, "K103"
label mut_V106 and name CA, "V106"
label mut_Y181 and name CA, "Y181"
label mut_Y188 and name CA, "Y188"
label mut_G190 and name CA, "G190"
label mut_P225 and name CA, "P225"
label mut_F227 and name CA, "F227"
label mut_M230 and name CA, "M230"
label mut_L234 and name CA, "L234"
label mut_Y318 and name CA, "Y318"
set label_size, 20
set label_color, black
set label_outline_color, white
set label_shadow_mode, 0
set label_bg_color, white
set label_bg_transparency, 0.25
set label_position, (0, 0, 3.0)

# ---- label nudges: edit these ----------------------------------------
# One line per residue. The three numbers are a screen-space offset in
# angstroms: x right, y up, z toward the viewer. The anchor sphere does
# not move, so a nudged label still marks its own residue.
#
# In the GUI you can also drag a label directly: enable editing with
# the Mouse menu (3-Button Editing), then ctrl-shift-drag the label.
# PyMOL stores the result in these same settings, so File > Save
# Session keeps whatever you do by hand.

# ---- scene -----------------------------------------------------------
bg_color white
set ray_opaque_background, 1
set cartoon_fancy_helices, 1
set antialias, 2
set ray_shadows, 0
deselect

# Suggested extras, left off so the first view stays readable:
#   show surface, NNIBP
#   set transparency, 0.5
#   show cartoon, p51   (p51 is already shown, semi-transparent)
set label_position, (0, 0, 3.0), mut_A98
set label_position, (0, 0, 3.0), mut_L100
set label_position, (0, 0, 3.0), mut_K103
set label_position, (0, 0, 3.0), mut_V106
set label_position, (0, 0, 3.0), mut_Y181
set label_position, (0, 0, 3.0), mut_Y188
set label_position, (0, 0, 3.0), mut_G190
set label_position, (0, 0, 3.0), mut_P225
set label_position, (0, 0, 3.0), mut_F227
set label_position, (0, 0, 3.0), mut_M230
set label_position, (0, 0, 3.0), mut_L234
set label_position, (0, 0, 3.0), mut_Y318

# Camera chosen by best_view(): of 4000 sampled directions, the
# one that maximises the smallest on-screen separation between
# labelled residues while keeping them clear of the ligand. It is
# a starting point, not a final answer -- P225/F227 and A98/Y318
# are only a few angstroms apart, so no direction separates them.
# Rotate as you like and re-save; nothing else depends on it.
set_view (\
  -0.416037, 0.909348, 0.000000,\
  -0.811820, -0.371417, 0.450552,\
  0.409709, 0.187446, 0.892750,\
  0.0, 0.0, -83.075,\
  9.303982, 11.424402, 16.513615,\
  70.488, 95.663, -20.0 )
