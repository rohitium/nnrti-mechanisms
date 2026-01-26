from __future__ import annotations

import os

from .require import require_module


def get_platform():
    openmm = require_module("openmm")
    platform_name = os.environ.get("OPENMM_PLATFORM", "").strip()
    if not platform_name:
        platform_name = (
            "Metal"
            if any(
                openmm.Platform.getPlatform(i).getName() == "Metal"
                for i in range(openmm.Platform.getNumPlatforms())
            )
            else "CPU"
        )
    platform = openmm.Platform.getPlatformByName(platform_name)
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
