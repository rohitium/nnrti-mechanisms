# TOC-left-panel style: faint NNIBP cartoon, pocket volume as a blue isomesh.
# No DOR. Run: /Applications/PyMOL.app/Contents/bin/pymol -cq <this file>
reinitialize
bg_color white
set ray_opaque_background, 0
set orthoscopic, 1
set ray_shadows, 0
set antialias, 2
set cartoon_transparency, 0.80
set mesh_width, 0.9
set stick_radius, 0.13
set mesh_quality, 2

python
from pymol import cmd
OUT = "results/analysis/apo_ensembles/pocket_render"
CASES = [("4NCG", "data/structures/apo_reference/4NCG.pdb", -3),
         ("1DLO", "data/structures/apo_reference/1DLO.pdb", -3),
         ("WT",   OUT + "/wt_apo_collapsed.pdb",             0)]
NNIBP = "97+98+100+103+104+105+176+178+185+186+187+224+226+231+315"

for tag, path, off in CASES:
    cmd.load(path, tag)
    if off:
        cmd.alter(tag, "resi = str(int(resi) + (%d))" % off)
        cmd.sort(tag)
    cmd.remove(f"{tag} and not polymer")            # drop DOR, waters, ions
    cmd.load(f"{OUT}/freevol_{tag}.dx", f"map_{tag}")
    chains = cmd.get_chains(tag)
    p66 = max((cmd.count_atoms(f"{tag} and chain {c} and name CA"), c) for c in chains)[1]
    cmd.select(f"pocket_{tag}", f"{tag} and chain {p66} and resi {NNIBP}")
    cmd.hide("everything", tag)
    cmd.show("cartoon", f"byres ({tag} within 14 of pocket_{tag})")
    # full residues, backbone included -- not just the side chain
    cmd.show("sticks", f"pocket_{tag} and not hydro")
    cmd.color("grey70", tag)
    cmd.color("palecyan", f"pocket_{tag}")
    cmd.util.cnc(f"pocket_{tag}")
    cmd.disable(tag)

# Superpose on the NNIBP, carrying each map with its own structure
for tag in ("1DLO", "WT"):
    cmd.align(f"pocket_{tag} and name CA", "pocket_4NCG and name CA")
    cmd.matrix_copy(tag, f"map_{tag}")

for tag, _, _ in CASES:
    cmd.isomesh(f"mesh_{tag}", f"map_{tag}", 0.5)
    cmd.color("marine", f"mesh_{tag}")
    cmd.disable(f"mesh_{tag}")

for tag, _, _ in CASES:
    cmd.enable(tag); cmd.enable(f"mesh_{tag}")
    cmd.orient("pocket_4NCG")
    cmd.zoom("pocket_4NCG", 7)
    cmd.ray(1100, 1100)
    cmd.png(f"{OUT}/pocket_{tag}.png", dpi=300)
    cmd.disable(tag); cmd.disable(f"mesh_{tag}")
    print(f"rendered {tag}")
python end
