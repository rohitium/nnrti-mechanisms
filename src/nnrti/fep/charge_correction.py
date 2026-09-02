#!/usr/bin/env python3
"""Analytical finite-size (net-charge) correction for charge-changing NEQ legs.

Background
----------
Our charge legs (wt_to_K103N: Lys+ -> Asn0; wt_to_G190E: Gly0 -> Glu-) change the
*system* net charge by delta_q = -1 during the alchemical switch. genion -neutral
neutralises the A-state, so the system goes q = 0 (state A) -> q = -1 (state B).
Under PME, a net-charged periodic cell is made neutral by a uniform background
("jellium"); the free energy of that periodic, background-neutralised system carries
a spurious box-size-dependent artifact relative to the true, macroscopic,
infinite-dilution result. Rocklin, Mobley, Dill & Hunenberger (J. Chem. Phys. 139,
184103, 2013; doi:10.1063/1.4826261) give the analytical correction.

Scope of THIS module (be honest about it)
------------------------------------------
It computes the **leading periodicity / net-charge self-energy term** (the dominant
analytical term; Hunenberger & McCammon, Biophys. Chem. 78, 69, 1999), which depends
only on the box and delta_q:

    dG_self(q) = xi_EW * q^2 / (8 pi eps0 eps_s L)

applied to the transformation as dG_cor = k * (q_B^2 - q_A^2) / L, with q_A = 0,
q_B = delta_q. xi_EW = -2.837297 is the cubic-lattice Wigner constant; we use the
cubic-equivalent edge L = V^(1/3) (the rhombic-dodecahedron xi differs by ~a few %,
far below the term's own size here).

It does NOT yet include the higher-order Rocklin terms: net-charge undersolvation
(~q^2/L^3, smaller and cancels similarly), the discrete-solvent / Galvani term
(q-linear, weakly density-dependent), and the residual integrated potential (RIP,
environment-dependent, the one term that does not cancel between holo and apo and
would need a Poisson-Boltzmann / explicit-potential pass). Those are documented as
residual uncertainty; see module-level NOTE in charge_leg_correction().

Why the leading term is essentially zero for OUR dG_bind
--------------------------------------------------------
dG_bind = dG_holo - dG_apo, and both phases undergo the SAME delta_q = -1 in nearly
identical boxes (holo V = 1722.6 nm^3, apo V = 1710.5 nm^3, ~0.7% apart). The per-leg
self-energy corrections (~0.04-0.05 kcal/mol each) therefore cancel to ~1e-4 kcal/mol
in the binding double difference. This is the expected result: relative charge-change
free energies in matched environments have tiny finite-size corrections. It also
vindicates dropping the co-alchemical ion, which injected ~20 kcal/mol of dissipation
to "fix" a <~0.3 kcal/mol artifact.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

# Physical constants (SI) + conversions
_E = 1.602176634e-19       # elementary charge, C
_EPS0 = 8.8541878128e-12   # vacuum permittivity, F/m
_NA = 6.02214076e23        # Avogadro, /mol
_KJ_PER_KCAL = 4.184
_XI_EW = -2.837297         # cubic-lattice (Wigner) self-energy constant

# TIP3P static dielectric ~= 97 (Hoefinger 2005 report ~94-97); experimental water is
# 78.4. The net dG_bind correction is insensitive to this (it cancels); the per-leg
# value scales as 1/eps_s. Default to the TIP3P model value we simulate with.
EPS_TIP3P = 97.0


def parse_gro_box(gro_path: str | Path) -> dict:
    """Return box volume (nm^3) and cubic-equivalent edge L (nm) from a .gro file.

    The last line of a GROMACS .gro is: v1x v2y v3z v1y v1z v2x v2z v3x v3y.
    Volume = det of the box matrix; for the lower-triangular GROMACS form that is
    v1x * v2y * v3z.
    """
    with open(gro_path) as handle:
        last = handle.readlines()[-1]
    f = [float(x) for x in last.split()]
    v1x, v2y, v3z = f[0], f[1], f[2]
    volume_nm3 = v1x * v2y * v3z
    l_eff_nm = volume_nm3 ** (1.0 / 3.0)
    return {"box_line": last.strip(), "volume_nm3": volume_nm3, "l_eff_nm": l_eff_nm}


def periodicity_self_energy(q_a: float, q_b: float, l_eff_nm: float, eps_s: float = EPS_TIP3P) -> float:
    """Leading net-charge periodicity correction, kcal/mol, to ADD to the raw dG.

    dG_cor = -xi_EW * (q_B^2 - q_A^2) * e^2 / (8 pi eps0 eps_s L)   [per charge]
    (the minus sign converts the periodic self-interaction into the correction to add
    to reach the infinite-system result). Charges q_a, q_b are in units of e.
    """
    l_m = l_eff_nm * 1e-9
    dq2 = (q_b ** 2) - (q_a ** 2)
    joules = -_XI_EW * dq2 * (_E ** 2) / (8.0 * math.pi * _EPS0 * eps_s * l_m)
    kcal_per_mol = joules * _NA / 1000.0 / _KJ_PER_KCAL
    return kcal_per_mol


def charge_leg_correction(
    holo_gro: str | Path,
    apo_gro: str | Path,
    *,
    delta_q: int = -1,
    eps_s: float = EPS_TIP3P,
) -> dict:
    """Per-phase leading corrections and the net effect on dG_bind (kcal/mol).

    NOTE: this is the leading periodicity term only. The non-cancelling residual
    (RIP / environment term) is NOT included here and is bounded separately; for a
    box this large with matched holo/apo cells it is expected to be a few tenths of
    a kcal/mol at most, i.e. within the per-leg statistical SEM.
    """
    holo = parse_gro_box(holo_gro)
    apo = parse_gro_box(apo_gro)
    dg_holo = periodicity_self_energy(0.0, delta_q, holo["l_eff_nm"], eps_s)
    dg_apo = periodicity_self_energy(0.0, delta_q, apo["l_eff_nm"], eps_s)
    return {
        "delta_q": delta_q,
        "eps_s": eps_s,
        "holo": {**holo, "dg_self_kcal": dg_holo},
        "apo": {**apo, "dg_self_kcal": dg_apo},
        "ddg_correction_kcal": dg_holo - dg_apo,  # add to raw dG_bind
    }


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Leading net-charge finite-size correction for charge legs.")
    p.add_argument("--holo-gro", required=True, type=Path)
    p.add_argument("--apo-gro", required=True, type=Path)
    p.add_argument("--delta-q", type=int, default=-1)
    p.add_argument("--eps-s", type=float, default=EPS_TIP3P)
    args = p.parse_args(argv)
    r = charge_leg_correction(args.holo_gro, args.apo_gro, delta_q=args.delta_q, eps_s=args.eps_s)
    print(f"delta_q={r['delta_q']}  eps_s={r['eps_s']}")
    for ph in ("holo", "apo"):
        d = r[ph]
        print(f"  {ph}: V={d['volume_nm3']:.1f} nm^3  L={d['l_eff_nm']:.3f} nm  "
              f"dG_self={d['dg_self_kcal']:+.4f} kcal/mol")
    print(f"  => net correction to dG_bind = {r['ddg_correction_kcal']:+.4f} kcal/mol")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
