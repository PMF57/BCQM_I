from __future__ import annotations
import json, sys, platform
from typing import Any, Dict

def save_metadata(path: str, meta: Dict[str, Any]) -> None:
    with open(path, "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, sort_keys=True)

def env_info() -> dict:
    return {"python": sys.version, "platform": platform.platform()}
