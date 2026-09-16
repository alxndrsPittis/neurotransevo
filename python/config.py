"""Path configuration shared by the Python scripts (mirrors R/config.R)."""
import os
from pathlib import Path

import yaml

ROOT = Path(os.environ.get("NEUROTRANSEVO_ROOT", Path(__file__).resolve().parents[1]))


def _merge(base, override):
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            _merge(base[key], value)
        else:
            base[key] = value
    return base


def load_config():
    cfg = yaml.safe_load((ROOT / "config.yml").read_text())
    local = ROOT / "config.local.yml"
    if local.exists():
        _merge(cfg, yaml.safe_load(local.read_text()) or {})
    paths = {}
    for key, value in cfg["paths"].items():
        value = os.environ.get(f"NEUROTRANSEVO_{key.upper()}", value)
        p = Path(value).expanduser()
        paths[key] = p if p.is_absolute() else ROOT / p
    cfg["paths"] = paths
    return cfg


CFG = load_config()


def results_path(*parts):
    p = CFG["paths"]["results"].joinpath(*parts)
    p.parent.mkdir(parents=True, exist_ok=True)
    return p
