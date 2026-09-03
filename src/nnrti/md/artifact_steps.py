from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path


def infer_json_steps(json_path: Path | None) -> int:
    if json_path is None or not json_path.exists():
        return 0
    try:
        payload = json.loads(json_path.read_text())
    except Exception:
        return 0
    value = payload.get("md_production_steps_completed", payload.get("md_production_steps", 0))
    try:
        return int(value or 0)
    except Exception:
        return 0


def infer_state_csv_steps(state_csv_path: Path | None) -> int:
    if state_csv_path is None or not state_csv_path.exists():
        return 0
    try:
        with state_csv_path.open(newline="") as handle:
            reader = csv.reader(handle)
            header = next(reader, None)
            if not header:
                return 0
            step_idx = None
            for idx, name in enumerate(header):
                if name in {'#"Step"', "Step"}:
                    step_idx = idx
                    break
            if step_idx is None:
                return 0
            max_step = 0
            for row in reader:
                if step_idx >= len(row):
                    continue
                raw = str(row[step_idx]).strip()
                if not raw:
                    continue
                try:
                    step = int(float(raw))
                except Exception:
                    continue
                if step > max_step:
                    max_step = step
            return max_step
    except Exception:
        return 0


def infer_state_csv_path(json_path: Path | None) -> Path | None:
    if json_path is None:
        return None
    if json_path.suffix.lower() != ".json":
        return None
    stem = json_path.stem.replace("_apo_rep", "_rep")
    return json_path.with_name(f"{stem}_md_state.csv")


@dataclass(frozen=True)
class StepInference:
    best_steps: int
    best_source: str
    json_steps: int
    state_csv_steps: int


def infer_best_steps(
    json_path: Path | None = None,
    state_csv_path: Path | None = None,
) -> StepInference:
    json_steps = infer_json_steps(json_path)
    if state_csv_path is None:
        state_csv_path = infer_state_csv_path(json_path)
    state_steps = infer_state_csv_steps(state_csv_path)

    if state_steps > json_steps:
        return StepInference(
            best_steps=state_steps,
            best_source="state_csv",
            json_steps=json_steps,
            state_csv_steps=state_steps,
        )
    if json_steps > 0:
        return StepInference(
            best_steps=json_steps,
            best_source="json",
            json_steps=json_steps,
            state_csv_steps=state_steps,
        )
    if state_steps > 0:
        return StepInference(
            best_steps=state_steps,
            best_source="state_csv",
            json_steps=json_steps,
            state_csv_steps=state_steps,
        )
    return StepInference(
        best_steps=0,
        best_source="missing",
        json_steps=json_steps,
        state_csv_steps=state_steps,
    )


@dataclass(frozen=True)
class ConsistencyStatus:
    json_path: str
    state_csv_path: str
    status: str
    json_steps: int
    state_csv_steps: int
    consistent: bool
    changed: bool


def read_json_payload(json_path: Path | None) -> dict:
    if json_path is None or not json_path.exists():
        return {}
    try:
        payload = json.loads(json_path.read_text())
    except Exception:
        return {}
    return payload if isinstance(payload, dict) else {}


def reconcile_json_with_state_csv(
    json_path: Path | None,
    state_csv_path: Path | None = None,
    *,
    write: bool = False,
    target_steps: int | None = None,
    allow_downgrade: bool = False,
) -> ConsistencyStatus:
    """Reconcile a run JSON against its state.csv.

    ``allow_downgrade`` guards the destructive direction. ``state.csv`` is often a
    **stale mid-slice dump** -- a resumed run can leave a file whose last Step is
    far below what the run actually completed -- so using it to LOWER
    ``md_production_steps_completed`` silently corrupts every time axis derived
    from that field (``result_collector`` builds frame timestamps from it).

    Observed 2026-08-17: an unconditional rewrite halved 25 run JSONs
    (e.g. G190A rep_01 50,000,000 -> 25,218,000) from state.csv files that the
    same tool had flagged ``state_csv_stale: true``. The analysis DCDs for those
    runs carry the canonical completed-100 ns fingerprint (180 fr / 34 MB), which
    the analysis DCD is authoritative over both the JSON
    and state.csv.

    Raising the count is safe and still automatic; lowering now requires an
    explicit opt-in.
    """
    payload = read_json_payload(json_path)
    if state_csv_path is None:
        state_csv_path = infer_state_csv_path(json_path)

    json_steps = infer_json_steps(json_path)
    state_csv_steps = infer_state_csv_steps(state_csv_path)
    status = str(payload.get("status", "")).lower()
    consistent = json_steps == state_csv_steps
    changed = False

    if (
        write
        and json_path is not None
        and json_path.exists()
        and state_csv_path is not None
        and state_csv_path.exists()
    ):
        updated = dict(payload)
        if target_steps is not None:
            updated["md_production_steps"] = int(target_steps)
        # Never LOWER the completion count from state.csv unless explicitly asked:
        # state.csv can be a stale mid-slice dump, and a downgrade silently
        # compresses every time axis derived from md_production_steps_completed.
        would_downgrade = 0 < state_csv_steps < json_steps
        if would_downgrade and not allow_downgrade:
            updated["state_csv_steps_observed"] = int(state_csv_steps)
            updated["state_csv_stale"] = True
        elif state_csv_steps > 0:
            updated["md_production_steps_completed"] = int(state_csv_steps)
            heating_steps = updated.get("md_heating_steps", 0)
            try:
                heating_steps_int = int(heating_steps or 0)
            except Exception:
                heating_steps_int = 0
            updated["md_total_steps"] = int(heating_steps_int + state_csv_steps)
        if updated != payload:
            json_path.write_text(json.dumps(updated, indent=2))
            payload = updated
            changed = True
            json_steps = infer_json_steps(json_path)
            status = str(payload.get("status", "")).lower()
            consistent = json_steps == state_csv_steps

    return ConsistencyStatus(
        json_path=str(json_path) if json_path is not None else "",
        state_csv_path=str(state_csv_path) if state_csv_path is not None else "",
        status=status,
        json_steps=json_steps,
        state_csv_steps=state_csv_steps,
        consistent=consistent,
        changed=changed,
    )
