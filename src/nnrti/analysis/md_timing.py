"""Infer production length for holo analysis DCDs.

JSON ``md_production_steps_completed`` and ``*_md_state.csv`` are often stale.
The worker writes ``status=ok`` after every SLURM slice, so a 12 h resume that
stops at 35 M steps freezes the metadata even if a later job finished the DCD.

DCD headers are also useless: OpenMM wrote ``nsavc=1``, ``DELTA=1.0 ps``.

What *is* reliable is the analysis DCD itself. Completed 100 ns runs cluster:

  180 frames, ~34.1 MB  — standard 100 ns product (Y188L, V106I, K103N, …)
  196–200 frames, ~37–38 MB
  360 frames, ~68.2 MB  — 2× frame density, still 100 ns (V106A+F227L, …)

A "50 ns" G190A DCD is the same size as a "100 ns" Y188L DCD. Treat those as
100 ns. Flag anything that does not match a known cluster.
"""
from __future__ import annotations

import json
import struct
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

DT_FS = 2.0
TARGET_NS = 100.0
TARGET_STEPS = 50_000_000

# n_frames -> typical size in MB for a completed 100 ns analysis DCD
_FINGERPRINTS_MB = {
    180: 34.12,
    196: 37.16,
    200: 37.91,
    360: 68.23,
}
_SIZE_TOL_MB = 0.35
_CANONICAL_FRAMES = 180
_CANONICAL_MB = 34.12
_MB_PER_FRAME_TOL = 0.003


@dataclass(frozen=True)
class TimingCall:
    mutation: str
    replicate: int
    n_frames: int | None
    dcd_mb: float | None
    json_ns: float | None
    state_csv_ns: float | None
    production_ns: float
    source: str
    note: str


def dcd_n_frames(dcd_path: Path) -> int | None:
    if not dcd_path.exists():
        return None
    with dcd_path.open("rb") as handle:
        header = handle.read(12)
    if len(header) < 12:
        return None
    if header[4:8] != b"CORD":
        return None
    return int(struct.unpack_from("<i", header, 8)[0])


def _json_ns(json_path: Path | None) -> float | None:
    if json_path is None or not json_path.exists():
        return None
    try:
        payload = json.loads(json_path.read_text())
    except Exception:
        return None
    steps = payload.get("md_production_steps_completed") or payload.get("md_production_steps")
    try:
        steps_i = int(steps or 0)
    except Exception:
        return None
    return steps_i * DT_FS / 1e6 if steps_i > 0 else None


def _state_csv_ns(csv_path: Path | None) -> float | None:
    if csv_path is None or not csv_path.exists():
        return None
    try:
        df = pd.read_csv(csv_path, usecols=lambda c: c in {'#"Step"', "Step"})
    except Exception:
        return None
    col = '#"Step"' if '#"Step"' in df.columns else "Step"
    steps = pd.to_numeric(df[col], errors="coerce").dropna()
    if steps.empty:
        return None
    return float(steps.max()) * DT_FS / 1e6


def _fingerprint_ns(n_frames: int | None, dcd_mb: float | None) -> str | None:
    if n_frames is None or dcd_mb is None:
        return None
    expected = _FINGERPRINTS_MB.get(int(n_frames))
    if expected is not None and abs(float(dcd_mb) - expected) <= _SIZE_TOL_MB:
        return "dcd_fingerprint_100ns"
    for frames, mb in _FINGERPRINTS_MB.items():
        if abs(float(dcd_mb) - mb) <= _SIZE_TOL_MB:
            return f"dcd_size_matches_{frames}fr_100ns"
    return None


def _short_same_stride_ns(n_frames: int | None, dcd_mb: float | None) -> float | None:
    """Fewer frames than the 180-frame 100 ns product, same bytes/frame → shorter run.

    G190E r3 is 139 frames / 26.4 MB vs siblings at 180 / 34.1 MB.
    Time = (n_frames-1) * (100 ns / 179).
    """
    if n_frames is None or dcd_mb is None or n_frames < 2 or n_frames >= _CANONICAL_FRAMES:
        return None
    expected_mb_per_frame = _CANONICAL_MB / float(_CANONICAL_FRAMES)
    actual = float(dcd_mb) / float(n_frames)
    if abs(actual - expected_mb_per_frame) > _MB_PER_FRAME_TOL:
        return None
    return float(n_frames - 1) * TARGET_NS / float(_CANONICAL_FRAMES - 1)


def infer_production_ns(
    *,
    dcd_path: Path | None,
    json_path: Path | None = None,
    state_csv_path: Path | None = None,
    mutation: str = "",
    replicate: int = 0,
) -> TimingCall:
    n_frames = dcd_n_frames(dcd_path) if dcd_path is not None else None
    dcd_mb = (dcd_path.stat().st_size / 1e6) if dcd_path is not None and dcd_path.exists() else None
    json_ns = _json_ns(json_path)
    csv_ns = _state_csv_ns(state_csv_path)
    fp = _fingerprint_ns(n_frames, dcd_mb)

    meta_vals = [v for v in (json_ns, csv_ns) if v is not None]
    meta_ns = max(meta_vals) if meta_vals else None

    if fp is not None:
        note = ""
        if meta_ns is not None and meta_ns < 95.0:
            note = f"metadata says {meta_ns:.1f} ns; DCD matches completed 100 ns cluster"
        return TimingCall(
            mutation=mutation,
            replicate=int(replicate),
            n_frames=n_frames,
            dcd_mb=dcd_mb,
            json_ns=json_ns,
            state_csv_ns=csv_ns,
            production_ns=TARGET_NS,
            source=fp,
            note=note,
        )

    short_ns = _short_same_stride_ns(n_frames, dcd_mb)
    if short_ns is not None:
        return TimingCall(
            mutation=mutation,
            replicate=int(replicate),
            n_frames=n_frames,
            dcd_mb=dcd_mb,
            json_ns=json_ns,
            state_csv_ns=csv_ns,
            production_ns=short_ns,
            source="dcd_short_same_stride",
            note=(
                f"{n_frames} frames / {dcd_mb:.1f} MB matches the 180-frame stride; "
                f"not a 100 ns product (metadata says {meta_ns:.1f} ns)"
                if meta_ns is not None
                else f"{n_frames} frames / {dcd_mb:.1f} MB matches the 180-frame stride; not a 100 ns product"
            ),
        )

    if meta_ns is not None and meta_ns >= 99.0:
        return TimingCall(
            mutation=mutation,
            replicate=int(replicate),
            n_frames=n_frames,
            dcd_mb=dcd_mb,
            json_ns=json_ns,
            state_csv_ns=csv_ns,
            production_ns=TARGET_NS,
            source="metadata_100ns",
            note="DCD frame count is non-canonical; metadata claims 100 ns",
        )

    # A 50 ns run at the 100 ns analysis stride is ~90–100 frames / ~17 MB.
    # Anything at or above the 180-frame / 34 MB cluster is not a short run.
    if n_frames is not None and dcd_mb is not None and n_frames >= 180 and dcd_mb >= 33.5:
        return TimingCall(
            mutation=mutation,
            replicate=int(replicate),
            n_frames=n_frames,
            dcd_mb=dcd_mb,
            json_ns=json_ns,
            state_csv_ns=csv_ns,
            production_ns=TARGET_NS,
            source="dcd_at_least_canonical_100ns",
            note=f"non-canonical {n_frames} frames / {dcd_mb:.1f} MB; too large to be a 50 ns DCD",
        )

    return TimingCall(
        mutation=mutation,
        replicate=int(replicate),
        n_frames=n_frames,
        dcd_mb=dcd_mb,
        json_ns=json_ns,
        state_csv_ns=csv_ns,
        production_ns=float(meta_ns) if meta_ns is not None else TARGET_NS,
        source="metadata_unverified",
        note="DCD does not match a known 100 ns fingerprint; check this run",
    )
