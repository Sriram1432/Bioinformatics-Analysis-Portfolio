"""Run Project 5 scripts in their actual dependency order."""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
# RNA modelling (03) must precede response follow-up (02), which optionally
# consumes held-out predictions produced by 03. Filenames are retained to avoid
# breaking existing references.
ORDER = (
    "01_prepare_data.py", "03_rna_models.py", "02_response_landscape.py",
    "04_multiomics_explainability.py", "05_cross_platform_validation.py",
    "06_nested_validation.py", "07_mapk_followup.py", "08_figures.py",
)
SCRIPTS = tuple(ROOT / name for name in ORDER)


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
