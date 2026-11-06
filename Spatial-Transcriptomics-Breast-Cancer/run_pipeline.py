"""Run the eight spatial-transcriptomics scripts in dependency order."""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = tuple(ROOT / f"{number:02d}_{name}.py" for number, name in enumerate((
    "prepare_data", "preprocess", "cluster", "characterise", "spatial_analysis",
    "spatial_genes", "pathways", "figures",
), start=1))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true", help="Print commands without executing analyses.")
    args = parser.parse_args()
    missing = [str(path) for path in SCRIPTS if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing pipeline scripts: {missing}")
    for script in SCRIPTS:
        command = [sys.executable, str(script)]
        print(" ".join(command), flush=True)
        if not args.dry_run:
            subprocess.run(command, cwd=ROOT, check=True)


if __name__ == "__main__":
    main()
