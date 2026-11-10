"""Load the pre-registration file ``config/analysis.yaml``."""

from __future__ import annotations

from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = REPO_ROOT / "config" / "analysis.yaml"


def load_config(path: str | Path | None = None) -> dict:
    """Read the YAML pre-registration (default: the repository's config/analysis.yaml)."""
    path = Path(path) if path is not None else DEFAULT_CONFIG
    with open(path) as fh:
        return yaml.safe_load(fh)
