#!/usr/bin/env python3
"""Plot WT-vs-mutation key DOR contact distances derived from 4NCG crystal contacts."""
from __future__ import annotations

import argparse
import json
import re
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _mutation_sort_key(m: str) -> tuple[int, str]:
    if m == "WT":
        return (0, m)
    if "+" in m:
        return (2, m)
    return (1, m)


def _remap_to_local_workspace(candidate: Path | None, repo_root: Path) -> Path | None:
    if candidate is None:
        return None
    if candidate.exists():
        return candidate
    marker = "nnrti-mechanisms/"
    text = str(candidate)
    if marker not in text:
        return candidate
    rel = text.split(marker, 1)[1]
    mapped = repo_root / rel
    if mapped.exists():
        return mapped
    return candidate


def _replicate_inputs(row: pd.Series, repo_root: Path) -> tuple[Path, Path]:
    data = json.loads(Path(row["output_json"]).read_text())
    topo = Path(str(data.get("analysis_topology_pdb") or "").strip())
    dcd = Path(str(data.get("analysis_dcd") or "").strip())
    topo = _remap_to_local_workspace(topo, repo_root)
    dcd = _remap_to_local_workspace(dcd, repo_root)
    if topo is None or dcd is None or not topo.exists() or not dcd.exists():
        raise FileNotFoundError(f"Missing analysis files for {row['mutation']} rep{int(row['replicate'])}")
    return topo, dcd


def _infer_total_ns_from_output_json(output_json_path: Path) -> float | None:
    state_csv = None
    m = re.match(r"^(.+)_rep(\d{2})\.json$", output_json_path.name)
    if m:
        safe = m.group(1)
        rep = int(m.group(2))
        state_csv = output_json_path.parent / f"{safe}_rep{rep:02d}_md_state.csv"
    if state_csv is None or not state_csv.exists():
        return None
    try:
        sdf = pd.read_csv(state_csv)
    except Exception:
        return None
    step_col = None
    for c in ('#"Step"', "Step"):
        if c in sdf.columns:
            step_col = c
            break
    if step_col is None or sdf.empty:
        return None
    steps = pd.to_numeric(sdf[step_col], errors="coerce").dropna()
    if steps.empty:
        return None
    return float(steps.max()) * 2.0 / 1_000_000.0


def _derive_contact_definitions(
    cif_path: Path,
    ligand_resname: str,
    ligand_auth_seq_id: int,
    protein_chain_id: str,
    polar_cutoff: float,
    hydrophobic_cutoff: float,
    top_polar: int,
    top_hydrophobic: int,
) -> pd.DataFrame:
    from src.utils.cif_parser import iter_cif_loops

    lines = cif_path.read_text().splitlines()
    atom_site_tags: list[str] | None = None
    atom_site_tokens: list[str] | None = None
    for tags, tokens in iter_cif_loops(lines):
        if "_atom_site.group_PDB" in tags and "_atom_site.auth_seq_id" in tags:
            atom_site_tags = tags
            atom_site_tokens = tokens
            break
    if atom_site_tags is None or atom_site_tokens is None:
        raise ValueError(f"Could not find _atom_site loop in {cif_path}")

    tags = atom_site_tags
    ncols = len(tags)
    idx = {t: tags.index(t) for t in tags}
    required = [
        "_atom_site.group_PDB",
        "_atom_site.type_symbol",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
        "_atom_site.auth_seq_id",
        "_atom_site.auth_comp_id",
        "_atom_site.auth_asym_id",
        "_atom_site.auth_atom_id",
    ]
    for r in required:
        if r not in idx:
            raise ValueError(f"Missing required mmCIF atom_site tag: {r}")

    lig_atoms: list[dict] = []
    prot_atoms: list[dict] = []
    for r0 in range(0, len(atom_site_tokens), ncols):
        row = atom_site_tokens[r0 : r0 + ncols]
        if len(row) < ncols:
            continue
        try:
            group = str(row[idx["_atom_site.group_PDB"]]).strip().upper()
            el = str(row[idx["_atom_site.type_symbol"]]).strip().upper()
            x = float(row[idx["_atom_site.Cartn_x"]])
            y = float(row[idx["_atom_site.Cartn_y"]])
            z = float(row[idx["_atom_site.Cartn_z"]])
            auth_seq = str(row[idx["_atom_site.auth_seq_id"]]).strip()
            comp = str(row[idx["_atom_site.auth_comp_id"]]).strip().upper()
            chain = str(row[idx["_atom_site.auth_asym_id"]]).strip()
            atom = str(row[idx["_atom_site.auth_atom_id"]]).strip()
        except Exception:
            continue

        if chain != str(protein_chain_id):
            continue
        if el == "H":
            continue

        if comp == str(ligand_resname).strip().upper() and auth_seq == str(int(ligand_auth_seq_id)):
            lig_atoms.append({"atom": atom, "el": el, "xyz": np.array([x, y, z], dtype=float)})
            continue

        if group == "ATOM":
            if not auth_seq.lstrip("+-").isdigit():
                continue
            resid = int(auth_seq)
            prot_atoms.append(
                {
                    "resid": resid,
                    "resname": comp,
                    "atom": atom,
                    "el": el,
                    "xyz": np.array([x, y, z], dtype=float),
                }
            )

    if not lig_atoms:
        raise ValueError(f"No heavy ligand atoms found for {ligand_resname} {protein_chain_id} {ligand_auth_seq_id} in {cif_path}")
    if not prot_atoms:
        raise ValueError(f"No protein atoms found in chain {protein_chain_id} in {cif_path}")

    polar_elements = {"N", "O", "S"}
    hydrophobic_lig_elements = {"C", "F", "CL", "BR", "I"}
    hydrophobic_prot_elements = {"C"}

    rows: list[dict] = []
    for la in lig_atoms:
        l_el = str(la["el"])
        l_xyz = la["xyz"]
        for pa in prot_atoms:
            p_el = str(pa["el"])
            p_xyz = pa["xyz"]
            dist = float(np.linalg.norm(l_xyz - p_xyz))

            category = None
            if l_el in polar_elements and p_el in polar_elements and dist <= float(polar_cutoff):
                category = "polar"
            elif l_el in hydrophobic_lig_elements and p_el in hydrophobic_prot_elements and dist <= float(hydrophobic_cutoff):
                category = "hydrophobic"
            if category is None:
                continue

            rows.append(
                {
                    "category": category,
                    "protein_chain": protein_chain_id,
                    "protein_resid_auth": int(pa["resid"]),
                    "protein_resname": str(pa["resname"]),
                    "protein_atom": str(pa["atom"]).strip(),
                    "ligand_resname": str(ligand_resname),
                    "ligand_atom": str(la["atom"]).strip(),
                    "distance_ref_angstrom": float(dist),
                }
            )

    if not rows:
        raise ValueError("No crystal contacts found with current thresholds.")
    df = pd.DataFrame(rows).sort_values("distance_ref_angstrom", ascending=True).reset_index(drop=True)

    # Keep shortest unique atom-pair contacts.
    df = df.drop_duplicates(
        subset=["category", "protein_resid_auth", "protein_atom", "ligand_atom"],
        keep="first",
    ).reset_index(drop=True)

    parts = []
    p = df[df["category"] == "polar"].head(max(0, int(top_polar)))
    h = df[df["category"] == "hydrophobic"].head(max(0, int(top_hydrophobic)))
    if not p.empty:
        parts.append(p)
    if not h.empty:
        parts.append(h)
    if not parts:
        raise ValueError("No contacts survived top-N filtering.")

    out = pd.concat(parts, ignore_index=True)
    out["contact_id"] = out.apply(
        lambda r: f"{r['category']}_{r['protein_resname']}{int(r['protein_resid_auth'])}_{r['protein_atom']}_{r['ligand_atom']}",
        axis=1,
    )
    out["contact_label"] = out.apply(
        lambda r: f"{r['category']}: {r['protein_resname']}{int(r['protein_resid_auth'])}:{r['protein_atom']} - {r['ligand_resname']}:{r['ligand_atom']}",
        axis=1,
    )
    return out


def _collect_system_rows(
    rows_df: pd.DataFrame,
    system_label: str,
    ligand_resname: str,
    frame_stride: int,
    repo_root: Path,
    contact_defs: pd.DataFrame,
    resid_offset: int,
) -> list[dict]:
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import distance_array

    out: list[dict] = []

    def _infer_element(atom_name: str) -> str:
        n = str(atom_name).strip().upper()
        if n.startswith("CL"):
            return "CL"
        if n.startswith("BR"):
            return "BR"
        if n.startswith("NA"):
            return "NA"
        if n.startswith("MG"):
            return "MG"
        if n.startswith("ZN"):
            return "ZN"
        return n[:1] if n else ""

    def _resolve_ligand_atom_group(u, ligand_resname: str, lig_atom_name: str, prot_atom_group):
        # Guard: if the protein selection is empty we cannot compute distances.
        if prot_atom_group.n_atoms == 0:
            return u.atoms[[]]

        # 1) Exact name
        ag = u.select_atoms(f"resname {ligand_resname} and name {lig_atom_name}")
        if ag.n_atoms == 1:
            return ag

        # 2) Prefix wildcard (handles x-suffixed renamed atoms like C1x/F1x)
        ag = u.select_atoms(f"resname {ligand_resname} and name {lig_atom_name}*")
        if ag.n_atoms == 1:
            return ag
        if ag.n_atoms > 1:
            # Pick nearest candidate to the protein atom at frame 0.
            d = distance_array(prot_atom_group.positions, ag.positions, box=u.dimensions).reshape(-1)
            if d.size == 0:
                return u.atoms[[]]
            return ag[[int(np.argmin(d))]]

        # 3) Element-based fallback for renamed atoms (e.g., 4NCG atom "F" -> trajectory "F1x")
        lig_all = u.select_atoms(f"resname {ligand_resname} and not name H*")
        if lig_all.n_atoms == 0:
            return lig_all
        target_el = _infer_element(lig_atom_name)
        cand_idx = []
        for i, a in enumerate(lig_all):
            el = str((getattr(a, "element", "") or a.name[:1]).upper())
            if el == target_el:
                cand_idx.append(i)
        if not cand_idx:
            return u.atoms[[]]
        cand = lig_all[cand_idx]
        d = distance_array(prot_atom_group.positions, cand.positions, box=u.dimensions).reshape(-1)
        if d.size == 0:
            return u.atoms[[]]
        return cand[[int(np.argmin(d))]]

    for _, row in rows_df.sort_values("replicate").iterrows():
        replicate = int(row["replicate"])
        try:
            topo, dcd = _replicate_inputs(row, repo_root)
            u = mda.Universe(str(topo), str(dcd))
        except Exception:
            continue

        # Resolve atom groups once per contact for this replicate.
        resolved = []
        for _, c in contact_defs.iterrows():
            resid = int(c["protein_resid_auth"]) + int(resid_offset)
            p_atom = str(c["protein_atom"])
            l_atom = str(c["ligand_atom"])
            contact_id = str(c["contact_id"])
            contact_label = str(c["contact_label"])
            category = str(c["category"])
            ref_d = float(c["distance_ref_angstrom"])

            p_sel = u.select_atoms(f"protein and resid {resid} and name {p_atom}")
            if p_sel.n_atoms == 0:
                continue
            l_sel = _resolve_ligand_atom_group(u, ligand_resname, l_atom, p_sel)
            if l_sel.n_atoms == 0:
                continue
            resolved.append((contact_id, contact_label, category, resid, p_atom, l_atom, ref_d, p_sel, l_sel))

        if not resolved:
            continue

        try:
            from MDAnalysis import transformations as trans

            prot = u.select_atoms("protein")
            anchor = prot if prot.n_atoms > 0 else u.atoms
            u.trajectory.add_transformations(
                trans.NoJump(check_continuity=False),
                trans.center_in_box(anchor, center="geometry", wrap=False),
            )
        except Exception:
            pass

        max_frame = max(1, len(u.trajectory) - 1)
        total_ns = _infer_total_ns_from_output_json(Path(str(row["output_json"])))
        if total_ns is None or not np.isfinite(total_ns) or total_ns <= 0:
            total_ns = 2.0

        for ts in u.trajectory[:: max(1, frame_stride)]:
            t_ns = (float(ts.frame) / float(max_frame)) * float(total_ns)
            for contact_id, contact_label, category, resid, p_atom, l_atom, ref_d, p_sel, l_sel in resolved:
                d = float(distance_array(p_sel.positions, l_sel.positions, box=u.dimensions).min())
                out.append(
                    {
                        "mutation": str(row["mutation"]),
                        "safe_label": str(row.get("safe_label", "")),
                        "system": system_label,
                        "replicate": replicate,
                        "time_ns": t_ns,
                        "contact_id": contact_id,
                        "contact_label": contact_label,
                        "category": category,
                        "protein_resid_traj": int(resid),
                        "protein_atom": p_atom,
                        "ligand_atom": l_atom,
                        "distance_ref_angstrom": ref_d,
                        "distance_angstrom": d,
                        "output_json": str(row.get("output_json", "")),
                        "analysis_dcd": str(dcd),
                        "analysis_topology_pdb": str(topo),
                    }
                )
    return out


def _interp_mean_trace(
    df: pd.DataFrame,
    x_col: str = "time_ns",
    y_col: str = "distance_angstrom",
    n_grid: int = 200,
):
    if df.empty:
        return None, None
    xmin = float(pd.to_numeric(df[x_col], errors="coerce").min())
    xmax = float(pd.to_numeric(df[x_col], errors="coerce").max())
    if not np.isfinite(xmin) or not np.isfinite(xmax) or xmax <= xmin:
        return None, None
    grid = np.linspace(xmin, xmax, n_grid)
    ys = []
    for _rep, grp in df.groupby("replicate"):
        g = grp.sort_values(x_col)
        x = g[x_col].to_numpy(dtype=float)
        y = g[y_col].to_numpy(dtype=float)
        if len(x) < 2:
            continue
        keep = np.r_[True, np.diff(x) > 0]
        x = x[keep]
        y = y[keep]
        if len(x) < 2:
            continue
        yi = np.interp(grid, x, y, left=np.nan, right=np.nan)
        yi[(grid < x.min()) | (grid > x.max())] = np.nan
        ys.append(yi)
    if not ys:
        return None, None
    return grid, np.nanmean(np.vstack(ys), axis=0)


def main() -> int:
    warnings.simplefilter("ignore", category=DeprecationWarning)

    parser = argparse.ArgumentParser(description="Plot crystal-derived DOR contact distances for all mutations vs WT.")
    parser.add_argument("--manifest", type=Path, default=Path("results/md_manifest.csv"))
    parser.add_argument("--cif", type=Path, default=Path("data/structures/4NCG.cif"))
    parser.add_argument("--ligand-resname", type=str, default="2KW")
    parser.add_argument("--ligand-auth-seq-id", type=int, default=601)
    parser.add_argument("--protein-chain", type=str, default="A")
    parser.add_argument("--frame-stride", type=int, default=1)
    parser.add_argument("--resid-offset", type=int, default=-3)
    parser.add_argument("--polar-cutoff", type=float, default=3.5)
    parser.add_argument("--hydrophobic-cutoff", type=float, default=4.5)
    parser.add_argument("--top-polar", type=int, default=4)
    parser.add_argument("--top-hydrophobic", type=int, default=8)
    parser.add_argument(
        "--mutations",
        type=str,
        default="",
        help="Optional comma-separated subset of mutations to process (e.g. V106I+F227C,V106A+L234I).",
    )
    parser.add_argument("--plots-dir", type=Path, default=Path("results/plots/dor_key_contacts"))
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=Path("results/dor_key_contacts_timeseries_all_mutations.csv"),
    )
    parser.add_argument(
        "--contact-defs-csv",
        type=Path,
        default=Path("results/dor_key_contact_definitions_4ncg.csv"),
    )
    args = parser.parse_args()

    if not args.manifest.exists():
        raise FileNotFoundError(args.manifest)
    if not args.cif.exists():
        raise FileNotFoundError(args.cif)

    repo_root = Path(__file__).resolve().parents[3]
    contact_defs = _derive_contact_definitions(
        cif_path=args.cif,
        ligand_resname=args.ligand_resname,
        ligand_auth_seq_id=args.ligand_auth_seq_id,
        protein_chain_id=args.protein_chain,
        polar_cutoff=float(args.polar_cutoff),
        hydrophobic_cutoff=float(args.hydrophobic_cutoff),
        top_polar=int(args.top_polar),
        top_hydrophobic=int(args.top_hydrophobic),
    )
    args.contact_defs_csv.parent.mkdir(parents=True, exist_ok=True)
    contact_defs.to_csv(args.contact_defs_csv, index=False)
    print(f"Wrote {args.contact_defs_csv} (n_contacts={len(contact_defs)})")

    mf = pd.read_csv(args.manifest)
    wt_df = mf[mf["mutation"] == "WT"].copy()
    mut_df = mf[mf["mutation"] != "WT"].copy()
    if wt_df.empty or mut_df.empty:
        raise ValueError("Manifest must contain WT and non-WT mutations.")

    muts = sorted(mut_df["mutation"].unique(), key=lambda m: _mutation_sort_key(str(m)))
    if str(args.mutations).strip():
        wanted = {m.strip() for m in str(args.mutations).split(",") if m.strip()}
        muts = [m for m in muts if str(m) in wanted]
        if not muts:
            raise ValueError("No requested mutations matched manifest entries.")
    all_rows: list[dict] = []

    mut_color = "#1f77b4"
    wt_color = "#444444"
    args.plots_dir.mkdir(parents=True, exist_ok=True)

    for mut in muts:
        mdf = mut_df[mut_df["mutation"] == mut].copy()
        wt_lookup = wt_df[wt_df["replicate"].isin(mdf["replicate"])].copy()
        if wt_lookup.empty:
            wt_lookup = wt_df.copy()

        rows = []
        rows.extend(
            _collect_system_rows(
                mdf,
                "Mutant",
                args.ligand_resname,
                args.frame_stride,
                repo_root,
                contact_defs,
                args.resid_offset,
            )
        )
        rows.extend(
            _collect_system_rows(
                wt_lookup,
                "WT",
                args.ligand_resname,
                args.frame_stride,
                repo_root,
                contact_defs,
                args.resid_offset,
            )
        )
        if not rows:
            continue

        block = pd.DataFrame(rows)
        all_rows.extend(block.to_dict(orient="records"))

        specs = contact_defs[["contact_id", "contact_label"]].drop_duplicates().values.tolist()
        nrows = len(specs)
        fig, axes = plt.subplots(nrows, 1, figsize=(10.0, max(3.2, 2.7 * nrows)), squeeze=False)
        axes_list = axes[:, 0].tolist()

        for i, (contact_id, contact_label) in enumerate(specs):
            ax = axes_list[i]
            sub = block[block["contact_id"] == str(contact_id)].copy()
            if sub.empty:
                ax.set_visible(False)
                continue

            for system, color in [("WT", wt_color), ("Mutant", mut_color)]:
                ss = sub[sub["system"] == system]
                for _rep, grp in ss.groupby("replicate"):
                    g = grp.sort_values("time_ns")
                    ax.plot(
                        g["time_ns"].to_numpy(dtype=float),
                        g["distance_angstrom"].to_numpy(dtype=float),
                        color=color,
                        alpha=0.25,
                        linewidth=0.8,
                    )

            for system, color in [("WT", wt_color), ("Mutant", mut_color)]:
                ss = sub[sub["system"] == system]
                x_mean, y_mean = _interp_mean_trace(ss)
                if x_mean is not None:
                    ax.plot(x_mean, y_mean, color=color, linewidth=2.0, alpha=0.95, label=system)
                if not ss.empty and ss["distance_angstrom"].notna().any():
                    mean_val = float(ss["distance_angstrom"].mean())
                    ax.axhline(mean_val, color=color, linestyle="--", linewidth=1.0, alpha=0.9)

            ref_d = float(pd.to_numeric(sub["distance_ref_angstrom"], errors="coerce").dropna().min())
            if np.isfinite(ref_d):
                ax.axhline(ref_d, color="#d62728", linestyle=":", linewidth=1.1, alpha=0.9, label="4NCG ref")

            xmax = float(pd.to_numeric(sub["time_ns"], errors="coerce").max())
            if np.isfinite(xmax) and xmax > 0:
                ax.set_xlim(0.0, xmax)
            ax.set_xlabel("Time (ns)")
            ax.set_ylabel("Distance (A)")
            ax.set_title(str(contact_label), fontsize=9)
            ax.grid(alpha=0.25, linestyle=":")
            ax.legend(frameon=False, fontsize=8, loc="best")

        fig.suptitle(f"{mut}: Crystal-Derived DOR Contact Distances (WT vs Mutant)", fontsize=11, fontweight="bold", y=0.995)
        fig.tight_layout()
        out_plot = args.plots_dir / f"{str(mut).replace('+', '_')}_dor_key_contacts_timeseries.png"
        fig.savefig(out_plot, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Wrote {out_plot}")

    if not all_rows:
        raise ValueError("No DOR-contact traces were generated.")

    out_df = pd.DataFrame(all_rows)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.output_csv, index=False)
    print(f"Wrote {args.output_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
