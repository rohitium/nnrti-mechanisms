#!/usr/bin/env python3
from pathlib import Path
import csv
from openmm import app, XmlSerializer

MANIFEST = Path("results/fep_manifest.csv")
TARGET_MUTATION = "V106A"
TARGET_REPLICATE = "1"

if not MANIFEST.exists():
    raise FileNotFoundError(f"Missing manifest: {MANIFEST}")

rows = []
with MANIFEST.open() as f:
    for r in csv.DictReader(f):
        if r["mutation"] == TARGET_MUTATION and r["replicate"] == TARGET_REPLICATE:
            rows.append(r)

if not rows:
    raise RuntimeError(f"No rows found for {TARGET_MUTATION} replicate {TARGET_REPLICATE}")

rows.sort(key=lambda r: r["leg"])

for r in rows:
    tid = r["task_id"]
    leg = r["leg"]
    pdb_txt = (r.get("prepared_topology_pdb") or "").strip()
    xml_txt = (r.get("prepared_system_xml") or "").strip()

    if not pdb_txt:
        raise RuntimeError(f"task {tid} ({leg}) has empty prepared_topology_pdb")
    if not xml_txt:
        raise RuntimeError(f"task {tid} ({leg}) has empty prepared_system_xml")

    pdb_path = Path(pdb_txt)
    xml_path = Path(xml_txt)

    if not pdb_path.exists():
        raise FileNotFoundError(f"task {tid} ({leg}) missing PDB: {pdb_path}")
    if not xml_path.exists():
        raise FileNotFoundError(f"task {tid} ({leg}) missing XML: {xml_path}")

    with pdb_path.open() as h:
        pdb = app.PDBFile(h)
    system = XmlSerializer.deserialize(xml_path.read_text())

    print(f"task {tid} {leg}: atoms={pdb.topology.getNumAtoms()} forces={system.getNumForces()} OK")

print("Preflight OK")
