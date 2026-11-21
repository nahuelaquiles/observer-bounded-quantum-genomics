"""Repository-relative paths shared by the analysis and figure scripts."""

import os
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = Path(os.environ.get("OBQG_DATA_DIR", REPO_ROOT / "data")).expanduser().resolve()
ASSET_DIR = Path(os.environ.get("OBQG_ASSET_DIR", REPO_ROOT / "assets")).expanduser().resolve()
RAW_DIR = Path(os.environ.get("OBQG_RAW_DIR", REPO_ROOT / "raw")).expanduser().resolve()
