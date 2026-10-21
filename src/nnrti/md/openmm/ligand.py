from __future__ import annotations

from pathlib import Path

from .require import require_module


#: Charge model used for the production MD and for every number in the main
#: text. Reviewer 2 (comment 5) asked how much this choice matters, so it is a
#: parameter rather than a constant -- but the default must stay "gasteiger",
#: because changing it silently would reparameterise the ligand away from the
#: MD that produced the trajectories being scored.
DEFAULT_CHARGE_METHOD = "gasteiger"


def _ensure_ambertools_registered() -> None:
    """Put AmberTools in the global toolkit registry if it is usable.

    openff-toolkit probes for antechamber and sqm at import time and silently
    omits the wrapper when they are not on PATH. Worker processes inherit a
    PATH without the conda environment's bin directory, so am1bcc resolves in
    an interactive shell and then fails inside the pool with "no registered
    toolkits can provide assign_partial_charges". Registering explicitly makes
    the behaviour the same in both places.
    """
    import os
    import sys
    from openff.toolkit.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY

    bindir = str(Path(sys.executable).parent)
    if bindir not in os.environ.get("PATH", "").split(os.pathsep):
        os.environ["PATH"] = bindir + os.pathsep + os.environ.get("PATH", "")
    try:
        from openff.toolkit.utils.ambertools_wrapper import AmberToolsToolkitWrapper
    except ImportError as exc:                                   # pragma: no cover
        raise RuntimeError("AmberTools wrapper unavailable; am1bcc needs it") from exc
    if any(isinstance(t, AmberToolsToolkitWrapper)
           for t in GLOBAL_TOOLKIT_REGISTRY.registered_toolkits):
        return
    wrapper = AmberToolsToolkitWrapper()
    if not wrapper.is_available():
        raise RuntimeError(
            "AmberTools reports unavailable; antechamber and sqm must be on PATH "
            f"(looked in {bindir})"
        )
    GLOBAL_TOOLKIT_REGISTRY.register_toolkit(wrapper)


def load_ligand_molecule(ligand_sdf: Path, charge_method: str = DEFAULT_CHARGE_METHOD):
    openff = require_module("openff.toolkit")
    mols = openff.topology.Molecule.from_file(
        str(ligand_sdf), allow_undefined_stereo=True
    )
    molecules = mols if isinstance(mols, list) else [mols]
    if not molecules:
        raise ValueError(f"No molecules found in {ligand_sdf}")
    if charge_method.startswith("am1"):
        _ensure_ambertools_registered()
    for mol in molecules:
        mol.assign_partial_charges(partial_charge_method=charge_method)
    return molecules[0]


#: SMIRNOFF release used to parameterise the ligand; verified against the
#: production system XMLs. See build_forcefield for why this is pinned.
SMIRNOFF_FORCEFIELD = "openff-2.0.0"


def build_forcefield(ligand_molecules,
                     charge_method: str = DEFAULT_CHARGE_METHOD) -> "openmm.app.ForceField":
    app = require_module("openmm.app")
    generators = require_module("openmmforcefields.generators")

    forcefield = app.ForceField(
        "amber14/protein.ff14SB.xml",
        "amber14/DNA.bsc1.xml",
        "amber14/tip3p.xml",
    )
    # Pinned explicitly. openmmforcefields <0.16 resolved forcefield=None to a
    # default release; 0.16 removed that fallback and raises instead, so the
    # version has to be named. openff-2.0.0 reproduces the ligand parameters
    # recorded in the production system XMLs bit-identically (sigma, epsilon and
    # charge all match to 0.0); openff-2.0.0 through 2.2.1 are equivalent for
    # doravirine, while 2.3.0 differs substantially. Changing this silently
    # reparameterises the ligand away from the MD that produced the trajectories.
    if charge_method.startswith("am1"):
        _ensure_ambertools_registered()
    generator = generators.SMIRNOFFTemplateGenerator(
        molecules=ligand_molecules,
        forcefield=SMIRNOFF_FORCEFIELD,
        template_generator_kwargs={"partial_charge_method": charge_method},
    )
    forcefield.registerTemplateGenerator(generator.generator)
    return forcefield
