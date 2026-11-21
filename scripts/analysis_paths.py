"""Repository-relative paths and portable loaders shared by the analysis and figure scripts."""

import os
from pathlib import Path
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = Path(os.environ.get("OBQG_DATA_DIR", REPO_ROOT / "data")).expanduser().resolve()
ASSET_DIR = Path(os.environ.get("OBQG_ASSET_DIR", REPO_ROOT / "assets")).expanduser().resolve()
RAW_DIR = Path(os.environ.get("OBQG_RAW_DIR", REPO_ROOT / "raw")).expanduser().resolve()

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def load_human_mane_promoters() -> pd.DataFrame:
    """Load the 4,570-promoter GRCh38 MANE Select cohort from OBQG_RAW_DIR if present,
    or reconstruct it identically from the portable gzipped sequence shards in data/."""
    feather_path = RAW_DIR / "human_grch38_mane_promoters.feather"
    if feather_path.exists():
        return pd.read_feather(feather_path)

    shards = [
        REPO_ROOT / "data" / "human_grch38_mane_promoters_part1.csv.gz",
        REPO_ROOT / "data" / "human_grch38_mane_promoters_part2.csv.gz",
    ]
    if all(p.exists() for p in shards):
        df = pd.concat([pd.read_csv(p) for p in shards], ignore_index=True)
        if "seq_sense_1001bp" not in df.columns:
            df["seq_sense_1001bp"] = [
                s if st == "+" else _revcomp(s)
                for s, st in zip(df["seq_plus_1001bp"], df["strand"])
            ]
        return df

    raise FileNotFoundError(
        f"Neither {feather_path} nor portable shards {shards[0]}..{shards[1]} were found."
    )
