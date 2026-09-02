"""Render GROMACS .mdp files for pmx NEQ runs."""

from __future__ import annotations

from pathlib import Path

from nnrti.fep.config import (
    NEQ_DT_PS,
    NEQ_EQUIL_NS,
    NEQ_TEMPERATURE_K,
    NEQ_WARMUP_PS,
    delta_lambda_for_switch,
    nsteps_for_time_ps,
)

_MDP_DIR = Path(__file__).resolve().parent / "mdp"


def _read_template(name: str) -> str:
    path = _MDP_DIR / name
    if not path.is_file():
        raise FileNotFoundError(path)
    return path.read_text()


def _write_mdp(path: Path, body: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(body.rstrip() + "\n")


def render_em_mdp(output: Path) -> None:
    _write_mdp(output, _read_template("em.mdp"))


def render_em_fep_mdp(*, output: Path, init_lambda: float) -> None:
    body = _read_template("em_fep.mdp")
    body = body.replace("@INIT_LAMBDA@", f"{init_lambda:.6f}")
    _write_mdp(output, body)


def render_npt_warmup_mdp(*, output: Path, init_lambda: float) -> None:
    nsteps = nsteps_for_time_ps(NEQ_WARMUP_PS)
    body = _read_template("npt_warmup.mdp")
    body = body.replace("@NSTEPS@", str(nsteps))
    body = body.replace("@INIT_LAMBDA@", f"{init_lambda:.6f}")
    body = body.replace("@REF_T@", f"{NEQ_TEMPERATURE_K:.2f}")
    _write_mdp(output, body)


def render_npt_eq_mdp(*, output: Path, init_lambda: float) -> None:
    nsteps = nsteps_for_time_ps(NEQ_EQUIL_NS * 1000.0)
    body = _read_template("npt_eq.mdp")
    body = body.replace("@NSTEPS@", str(nsteps))
    body = body.replace("@INIT_LAMBDA@", f"{init_lambda:.6f}")
    body = body.replace("@REF_T@", f"{NEQ_TEMPERATURE_K:.2f}")
    _write_mdp(output, body)


def render_nonequil_mdp(
    *,
    output: Path,
    init_lambda: float,
    switch_ps: float,
) -> None:
    nsteps = nsteps_for_time_ps(switch_ps)
    delta = delta_lambda_for_switch(switch_ps)
    if init_lambda >= 0.5:
        delta = -delta
    body = _read_template("nonequil.mdp")
    body = body.replace("@NSTEPS@", str(nsteps))
    body = body.replace("@INIT_LAMBDA@", f"{init_lambda:.6f}")
    body = body.replace("@DELTA_LAMBDA@", f"{delta:.8e}")
    body = body.replace("@REF_T@", f"{NEQ_TEMPERATURE_K:.2f}")
    _write_mdp(output, body)
