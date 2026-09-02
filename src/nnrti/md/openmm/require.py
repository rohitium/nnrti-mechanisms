from __future__ import annotations


def require_module(module_name: str):
    try:
        import importlib

        return importlib.import_module(module_name)
    except Exception as exc:  # pragma: no cover
        raise RuntimeError(
            f"Missing dependency '{module_name}'. Install required packages and retry."
        ) from exc
