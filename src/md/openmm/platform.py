from __future__ import annotations

import os
import logging

from .require import require_module


logger = logging.getLogger(__name__)


def get_platform():
    openmm = require_module("openmm")
    available = [
        openmm.Platform.getPlatform(i).getName()
        for i in range(openmm.Platform.getNumPlatforms())
    ]
    requested = os.environ.get("OPENMM_PLATFORM", "").strip()

    cuda_visible = os.environ.get("CUDA_VISIBLE_DEVICES")
    logger.info(
        "Platform selection: available=%s, requested=%s, CUDA_VISIBLE_DEVICES=%s",
        available, requested or "(auto)", cuda_visible,
    )

    # Preference order if not explicitly requested.
    if requested:
        candidates = [requested] + [p for p in ("CUDA", "OpenCL", "Metal", "CPU") if p != requested]
    else:
        candidates = ["CUDA", "OpenCL", "Metal", "CPU"]

    platform = None
    platform_name = None
    last_exc = None
    for name in candidates:
        if name not in available:
            continue
        try:
            platform = openmm.Platform.getPlatformByName(name)
            platform_name = name
            break
        except Exception as exc:  # pragma: no cover - defensive
            last_exc = exc
            continue

    if platform is None:
        if last_exc is not None:
            raise last_exc
        raise RuntimeError(f"No compatible OpenMM platform found. Available: {available}")

    if requested and platform_name != requested:
        logging.warning(
            "Requested OPENMM_PLATFORM=%s not usable; falling back to %s",
            requested,
            platform_name,
        )

    properties = {}
    if platform_name == "CPU":
        threads = os.environ.get("OPENMM_CPU_THREADS")
        if threads:
            properties["Threads"] = threads
    if platform_name in {"OpenCL", "Metal", "CUDA"}:
        device_index = os.environ.get("OPENMM_DEVICE_INDEX")
        if device_index:
            properties["DeviceIndex"] = device_index
    return platform, properties
