"""Runtime patches for Perses behavior we rely on in hybrid mutation prep."""

from __future__ import annotations


def patch_perses_proline_support() -> None:
    """Allow PRO in Perses charge-difference handling.

    Perses deliberately omits PRO from ``PolymerProposalEngine._aminos`` because
    of an OpenMM template issue, but ``PointMutationEngine._handle_charge_changes``
    still asserts both residue names are members before computing a (usually zero)
    charge difference. PRO->X mutations such as P225H therefore fail unless PRO
    is treated as a neutral amino acid here.
    """
    from perses.rjmc.topology_proposal import PolymerProposalEngine

    if "PRO" in PolymerProposalEngine._aminos:
        return
    aminos = list(PolymerProposalEngine._aminos)
    aminos.append("PRO")
    PolymerProposalEngine._aminos = aminos
    PolymerProposalEngine._neutral_aminos = PolymerProposalEngine._get_neutrals(
        PolymerProposalEngine._aminos,
        PolymerProposalEngine._positive_aminos,
        PolymerProposalEngine._negative_aminos,
    )
