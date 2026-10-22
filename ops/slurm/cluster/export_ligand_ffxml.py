#!/usr/bin/env python3
"""Write DOR's OpenMM force-field XML so the cluster never needs openff-toolkit.

openff-toolkit is distributed through conda-forge only; it is not on PyPI, and
Sherlock has no conda. The dependency exists solely to turn DOR into OpenMM
parameters, and DOR is the same molecule in every system, so the parameters can
be generated once here and shipped as an ffxml that plain OpenMM can load.

Usage
-----
    PYTHONPATH=src python ops/slurm/cluster/export_ligand_ffxml.py \
        data/ligands/dor_am1bcc.sdf data/ligands/dor_am1bcc.xml
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))


def main() -> int:
    from openmmforcefields.generators import SMIRNOFFTemplateGenerator
    from nnrti.md.openmm.ligand import SMIRNOFF_FORCEFIELD, load_ligand_molecule

    src, dst = Path(sys.argv[1]), Path(sys.argv[2])
    charge_method = sys.argv[3] if len(sys.argv) > 3 else None
    if charge_method:
        mol = load_ligand_molecule(src, charge_method)
    else:
        # The SDF already carries its charges. load_ligand_molecule would call
        # assign_partial_charges and overwrite them, which silently produced
        # two identical ffxmls the first time this was run.
        from openff.toolkit.topology import Molecule
        mol = Molecule.from_file(str(src), allow_undefined_stereo=True)
        if mol.partial_charges is None:
            raise SystemExit(f"{src} carries no partial charges; "
                             "pass a charge method as the third argument")
    import numpy as np
    q = np.array([float(x.m) for x in mol.partial_charges])
    print(f"ligand charges: mean |q| {np.abs(q).mean():.4f} e, "
          f"net {q.sum():+.2e}")
    gen = SMIRNOFFTemplateGenerator(molecules=[mol],
                                    forcefield=SMIRNOFF_FORCEFIELD)
    ffxml = gen.generate_residue_template(mol)
    if ffxml is None:
        ffxml = gen._ffxml_cache.get(mol.to_smiles())
    if not ffxml:
        raise SystemExit("generator produced no ffxml")
    dst.write_text(ffxml)

    # load it back through plain OpenMM, which is all the cluster will have
    import openmm.app as app
    ff = app.ForceField("amber14/protein.ff14SB.xml", "amber14/tip3p.xml", str(dst))
    print(f"wrote {dst} ({len(ffxml)} chars); "
          f"plain OpenMM ForceField loads it, {len(ff._templates)} templates")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
