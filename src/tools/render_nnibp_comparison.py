#!/usr/bin/env python3
"""
Render side-by-side NNIBP comparison using NGL.js directly.

Outputs an HTML file with interactive 3D view.

Usage:
    arch -arm64 uv run python -m src.tools.render_nnibp_comparison
"""
from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align


def main():
    ROOT = Path(__file__).resolve().parents[2]

    # Paths
    wt_pdb = ROOT / "data" / "generated" / "rpv" / "wt" / "rep_01" / "wt_minimized_rep01.pdb"
    mut_pdb = ROOT / "data" / "generated" / "rpv" / "Y181C_K101E" / "rep_01" / "mut_minimized_Y181C_K101E_rep01.pdb"
    output_html = ROOT / "results" / "plots" / "nnibp_comparison_y181c_k101e.html"

    # Config
    LIGAND_RESNAME = "T27"
    DRM_RESIS = [94, 174]  # K101, Y181 in PDB numbering
    NNIBP_RESIS = [93, 94, 96, 99, 101, 172, 174, 176, 181, 183, 220, 222, 223]
    POCKET_RADIUS = 8.0
    SHIFT_X = 25.0

    # Metrics
    WT_VOL, MUT_VOL = 1151.6, 1205.5
    WT_CONTACTS, MUT_CONTACTS = 552, 506

    # Colors
    WT_COLOR = "0x4CAF50"
    MUT_COLOR = "0xFF5722"
    DRM_COLOR = "0xE91E63"
    LIG_COLOR = "0x2196F3"

    print("Loading structures...")
    u_wt = mda.Universe(str(wt_pdb))
    u_mut = mda.Universe(str(mut_pdb))

    # Align
    align.alignto(u_mut, u_wt, select="protein and backbone")

    # Build pocket selections (chain A only - p66 subunit contains the NNIBP)
    lig_sel = f"resname {LIGAND_RESNAME}"
    nnibp_sel = "resid " + " ".join(map(str, NNIBP_RESIS))
    pocket_sel = f"(chainID A and protein and (around {POCKET_RADIUS} ({lig_sel}))) or ({lig_sel}) or (chainID A and protein and ({nnibp_sel}))"

    pocket_wt = u_wt.select_atoms(pocket_sel)
    pocket_mut = u_mut.select_atoms(pocket_sel)

    # Write pocket PDBs
    tmp_dir = Path(tempfile.mkdtemp(prefix="nnrti_"))
    out_wt = tmp_dir / "wt_pocket.pdb"
    out_mut = tmp_dir / "mut_pocket.pdb"
    out_mut_shift = tmp_dir / "mut_pocket_shifted.pdb"

    pocket_wt.write(str(out_wt))
    pocket_mut.write(str(out_mut))

    # Shift MUT
    u_mut_pocket = mda.Universe(str(out_mut))
    u_mut_pocket.atoms.positions += np.array([SHIFT_X, 0.0, 0.0])
    u_mut_pocket.atoms.write(str(out_mut_shift))

    # Calculate centroid distances for labels
    # WT distances
    wt_pocket = mda.Universe(str(out_wt))
    wt_lig = wt_pocket.select_atoms(f"resname {LIGAND_RESNAME}")
    wt_k101 = wt_pocket.select_atoms("resid 94")
    wt_y181 = wt_pocket.select_atoms("resid 174")
    wt_lig_center = wt_lig.center_of_geometry() if wt_lig.n_atoms > 0 else np.zeros(3)
    wt_k101_center = wt_k101.center_of_geometry() if wt_k101.n_atoms > 0 else np.zeros(3)
    wt_y181_center = wt_y181.center_of_geometry() if wt_y181.n_atoms > 0 else np.zeros(3)

    wt_dist_k101 = np.linalg.norm(wt_k101_center - wt_lig_center) if wt_k101.n_atoms > 0 and wt_lig.n_atoms > 0 else 0
    wt_dist_y181 = np.linalg.norm(wt_y181_center - wt_lig_center) if wt_y181.n_atoms > 0 and wt_lig.n_atoms > 0 else 0

    # MUT distances (use unshifted for accurate measurement)
    mut_pocket = mda.Universe(str(out_mut))
    mut_lig = mut_pocket.select_atoms(f"resname {LIGAND_RESNAME}")
    mut_k101 = mut_pocket.select_atoms("resid 94")
    mut_y181 = mut_pocket.select_atoms("resid 174")
    mut_lig_center = mut_lig.center_of_geometry() if mut_lig.n_atoms > 0 else np.zeros(3)
    mut_k101_center = mut_k101.center_of_geometry() if mut_k101.n_atoms > 0 else np.zeros(3)
    mut_y181_center = mut_y181.center_of_geometry() if mut_y181.n_atoms > 0 else np.zeros(3)

    mut_dist_k101 = np.linalg.norm(mut_k101_center - mut_lig_center) if mut_k101.n_atoms > 0 and mut_lig.n_atoms > 0 else 0
    mut_dist_y181 = np.linalg.norm(mut_y181_center - mut_lig_center) if mut_y181.n_atoms > 0 and mut_lig.n_atoms > 0 else 0

    # Pocket width (distance between K101 and Y181)
    wt_k101_ca = wt_pocket.select_atoms("name CA and resid 94")
    wt_y181_ca = wt_pocket.select_atoms("name CA and resid 174")
    mut_k101_ca = mut_pocket.select_atoms("name CA and resid 94")
    mut_y181_ca = mut_pocket.select_atoms("name CA and resid 174")
    wt_pocket_width = np.linalg.norm(wt_k101_ca.positions[0] - wt_y181_ca.positions[0]) if wt_k101_ca.n_atoms > 0 and wt_y181_ca.n_atoms > 0 else 0
    mut_pocket_width = np.linalg.norm(mut_k101_ca.positions[0] - mut_y181_ca.positions[0]) if mut_k101_ca.n_atoms > 0 and mut_y181_ca.n_atoms > 0 else 0

    print(f"WT: K101-RPV={wt_dist_k101:.1f}Å, Y181-RPV={wt_dist_y181:.1f}Å, pocket={wt_pocket_width:.1f}Å")
    print(f"MUT: E101-RPV={mut_dist_k101:.1f}Å, C181-RPV={mut_dist_y181:.1f}Å, pocket={mut_pocket_width:.1f}Å")

    def _append_centroid_atoms(pdb_content: str, centroids: list[tuple[str, np.ndarray]]) -> str:
        lines = pdb_content.rstrip("\n").splitlines()
        max_serial = max(int(line[6:11]) for line in lines if line.startswith(("ATOM", "HETATM")))
        insert_idx = len(lines) - 1 if lines and lines[-1].startswith("END") else len(lines)
        serial = max_serial + 1
        for name, pos in centroids:
            lines.insert(
                insert_idx,
                f"HETATM{serial:5d} {name:>4s} CEN A 999    {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}  1.00  0.00           C  ",
            )
            insert_idx += 1
            serial += 1
        return "\n".join(lines) + "\n"

    # Read PDB content and append centroid atoms
    with open(out_wt) as f:
        wt_pdb_raw = f.read()
    wt_pdb_raw = _append_centroid_atoms(
        wt_pdb_raw,
        [
            ("LIG", wt_lig_center),
            ("K101", wt_k101_center),
            ("Y181", wt_y181_center),
        ],
    )
    with open(out_mut_shift) as f:
        mut_pdb_raw = f.read()
    mut_pdb_raw = _append_centroid_atoms(
        mut_pdb_raw,
        [
            ("LIG", mut_lig_center + np.array([SHIFT_X, 0.0, 0.0])),
            ("K101", mut_k101_center + np.array([SHIFT_X, 0.0, 0.0])),
            ("Y181", mut_y181_center + np.array([SHIFT_X, 0.0, 0.0])),
        ],
    )
    wt_pdb_content = wt_pdb_raw.replace("\n", "\\n").replace("'", "\\'")
    mut_pdb_content = mut_pdb_raw.replace("\n", "\\n").replace("'", "\\'")

    # NGL selections (chain A only)
    nnibp_ngl = " or ".join([f"{r}:A" for r in NNIBP_RESIS])
    drm_ngl = " or ".join([f"{r}:A" for r in DRM_RESIS])

    print("Building HTML...")
    output_html.parent.mkdir(parents=True, exist_ok=True)

    html = f"""<!DOCTYPE html>
<html>
<head>
    <title>NNIBP Comparison: WT vs Y181C+K101E</title>
    <script src="https://unpkg.com/ngl@2.0.0-dev.37/dist/ngl.js"></script>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        html, body {{ width: 100%; height: 100%; overflow: hidden; font-family: Arial, sans-serif; }}
        #viewport {{ width: 100%; height: 100%; }}
        #metrics {{
            position: fixed; top: 10px; left: 50%; transform: translateX(-50%); z-index: 100;
            display: flex; gap: 20px;
        }}
        .metric-box {{
            padding: 10px 15px; border-radius: 8px; background: rgba(255,255,255,0.95);
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        .wt {{ border: 2px solid #4CAF50; }}
        .mut {{ border: 2px solid #FF5722; }}
        .wt h3 {{ color: #4CAF50; margin: 0 0 5px 0; font-size: 14px; }}
        .mut h3 {{ color: #FF5722; margin: 0 0 5px 0; font-size: 14px; }}
        .metric-box div {{ font-size: 13px; }}
        .delta {{ color: red; }}
        #legend {{
            position: fixed; bottom: 10px; left: 50%; transform: translateX(-50%); z-index: 100;
            background: rgba(255,255,255,0.95); padding: 8px 15px; border-radius: 5px;
            font-size: 14px; box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        #legend span {{ margin: 0 10px; }}
    </style>
</head>
<body>
    <div id="metrics">
        <div class="metric-box wt">
            <h3>Wild Type</h3>
            <div>Vol: <b>{WT_VOL:.1f} Å³</b> | Contacts: <b>{WT_CONTACTS}</b></div>
        </div>
        <div class="metric-box mut">
            <h3>Y181C + K101E</h3>
            <div>Vol: <b>{MUT_VOL:.1f} Å³</b> <span class="delta">(Δ +{MUT_VOL - WT_VOL:.1f})</span> | Contacts: <b>{MUT_CONTACTS}</b> <span class="delta">(Δ {MUT_CONTACTS - WT_CONTACTS:+d})</span></div>
        </div>
    </div>
    <div id="legend">
        <span style="color: #4CAF50;">■</span> WT
        <span style="color: #FF5722;">■</span> MUT
        <span style="color: #2196F3;">■</span> RPV
        <span style="color: #E91E63;">■</span> DRM sites
    </div>
    <div id="viewport"></div>
    <script>
        document.addEventListener("DOMContentLoaded", function() {{
            var stage = new NGL.Stage("viewport", {{ backgroundColor: "white" }});

            window.addEventListener("resize", function() {{ stage.handleResize(); }});

            // WT structure
            var wtPdb = '{wt_pdb_content}';
            var wtBlob = new Blob([wtPdb], {{ type: 'text/plain' }});

            stage.loadFile(wtBlob, {{ ext: "pdb", name: "wt" }}).then(function(comp) {{
                // Cartoon for protein backbone
                comp.addRepresentation("cartoon", {{
                    sele: "protein",
                    color: {WT_COLOR},
                    opacity: 0.5
                }});
                // DRM residues as ball+stick
                comp.addRepresentation("ball+stick", {{
                    sele: "94 or 174",
                    color: {DRM_COLOR},
                    multipleBond: true
                }});
                // Ligand as ball+stick
                comp.addRepresentation("ball+stick", {{
                    sele: "T27",
                    color: {LIG_COLOR},
                    multipleBond: true
                }});
                // DRM Labels using format
                comp.addRepresentation("label", {{
                    sele: "94 and .CA",
                    labelType: "format",
                    labelFormat: "K101",
                    color: "white",
                    labelSize: 2.5,
                    showBackground: true,
                    backgroundColor: "#E91E63",
                    backgroundOpacity: 0.9,
                    backgroundMargin: 3,
                    attachment: "middle-center"
                }});
                comp.addRepresentation("label", {{
                    sele: "174 and .CA",
                    labelType: "format",
                    labelFormat: "Y181",
                    color: "white",
                    labelSize: 2.5,
                    showBackground: true,
                    backgroundColor: "#E91E63",
                    backgroundOpacity: 0.9,
                    backgroundMargin: 3,
                    attachment: "middle-center"
                }});
                // RPV label
                comp.addRepresentation("label", {{
                    sele: "T27 and .C1x",
                    labelType: "format",
                    labelFormat: "RPV",
                    color: "white",
                    labelSize: 2.5,
                    showBackground: true,
                    backgroundColor: "#2196F3",
                    backgroundOpacity: 0.9,
                    backgroundMargin: 3,
                    attachment: "middle-center"
                }});
                // Distance lines between centroids
                comp.addRepresentation("distance", {{
                    atomPair: [
                        ["CEN and .K101", "CEN and .LIG"],
                        ["CEN and .Y181", "CEN and .LIG"]
                    ],
                    color: "black",
                    labelSize: 3.0,
                    labelColor: "black",
                    labelBackground: true,
                    labelBackgroundColor: "white",
                    labelBackgroundOpacity: 0.9,
                    labelUnit: "angstrom"
                }});
            }});

            // MUT structure (shifted)
            var mutPdb = '{mut_pdb_content}';
            var mutBlob = new Blob([mutPdb], {{ type: 'text/plain' }});

            stage.loadFile(mutBlob, {{ ext: "pdb", name: "mut" }}).then(function(comp) {{
                // Cartoon for protein backbone
                comp.addRepresentation("cartoon", {{
                    sele: "protein",
                    color: {MUT_COLOR},
                    opacity: 0.5
                }});
                // DRM residues as ball+stick
                comp.addRepresentation("ball+stick", {{
                    sele: "94 or 174",
                    color: {DRM_COLOR},
                    multipleBond: true
                }});
                // Ligand as ball+stick
                comp.addRepresentation("ball+stick", {{
                    sele: "T27",
                    color: {LIG_COLOR},
                    multipleBond: true
                }});
                // DRM Labels using format
                comp.addRepresentation("label", {{
                    sele: "94 and .CA",
                    labelType: "format",
                    labelFormat: "E101",
                    color: "white",
                    labelSize: 2.5,
                    showBackground: true,
                    backgroundColor: "#E91E63",
                    backgroundOpacity: 0.9,
                    backgroundMargin: 3,
                    attachment: "middle-center"
                }});
                comp.addRepresentation("label", {{
                    sele: "174 and .CA",
                    labelType: "format",
                    labelFormat: "C181",
                    color: "white",
                    labelSize: 2.5,
                    showBackground: true,
                    backgroundColor: "#E91E63",
                    backgroundOpacity: 0.9,
                    backgroundMargin: 3,
                    attachment: "middle-center"
                }});
                // RPV label
                comp.addRepresentation("label", {{
                    sele: "T27 and .C1x",
                    labelType: "format",
                    labelFormat: "RPV",
                    color: "white",
                    labelSize: 2.5,
                    showBackground: true,
                    backgroundColor: "#2196F3",
                    backgroundOpacity: 0.9,
                    backgroundMargin: 3,
                    attachment: "middle-center"
                }});
                // Distance lines between centroids
                comp.addRepresentation("distance", {{
                    atomPair: [
                        ["CEN and .K101", "CEN and .LIG"],
                        ["CEN and .Y181", "CEN and .LIG"]
                    ],
                    color: "black",
                    labelSize: 3.0,
                    labelColor: "black",
                    labelBackground: true,
                    labelBackgroundColor: "white",
                    labelBackgroundOpacity: 0.9,
                    labelUnit: "angstrom"
                }});

                // Center view after both loaded
                stage.autoView();
            }});
        }});
    </script>
</body>
</html>
"""

    with open(output_html, "w") as f:
        f.write(html)

    print(f"Done! Open in browser: {output_html}")


if __name__ == "__main__":
    main()
