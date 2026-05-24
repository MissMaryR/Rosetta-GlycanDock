"""
GlycanDock Prepack Submission Script

Finds all [GH]/[LIG]/1prepack/ directories under the given GlycanDock root,
skips any that already have results/packed.pdb, and submits the rest via
the submit.sh already present in each 1prepack/ folder.

Usage:
    python submit_all_prepack.py               # uses current working directory
    python submit_all_prepack.py <directory>   # uses specified GlycanDock root
"""

import sys
import subprocess
from pathlib import Path
import argparse


def main():
    parser = argparse.ArgumentParser(
        description="Submit all GlycanDock prepack jobs (Stage 0)"
    )
    parser.add_argument(
        "directory",
        nargs="?",
        default=".",
        help="GlycanDock root directory (default: current directory)",
    )
    args = parser.parse_args()

    root = Path(args.directory).resolve()

    if not root.exists():
        print(f"Error: directory {root} does not exist")
        sys.exit(1)
    if not root.is_dir():
        print(f"Error: {root} is not a directory")
        sys.exit(1)

    # Discover all 1prepack/ directories (pattern: GH/LIG/1prepack/)
    prepack_dirs = sorted(root.glob("*/*/1prepack"))

    if not prepack_dirs:
        print(f"No [GH]/[LIG]/1prepack/ directories found under {root}")
        sys.exit(1)

    print(f"GlycanDock root : {root}")
    print(f"Prepack folders : {len(prepack_dirs)} found")
    print()

    submitted = []
    skipped   = []
    failed    = []

    for prepack_dir in prepack_dirs:
        lig = prepack_dir.parent.name
        gh  = prepack_dir.parent.parent.name
        label = f"{gh}/{lig}"

        submit_sh = prepack_dir / "submit.sh"
        packed    = prepack_dir / "results" / "packed.pdb"

        if packed.exists():
            print(f"  SKIP   {label}  (results/packed.pdb already exists)")
            skipped.append(label)
            continue

        if not submit_sh.exists():
            print(f"  ERROR  {label}  (no submit.sh found in {prepack_dir})")
            failed.append(label)
            continue

        try:
            result = subprocess.run(
                ["sbatch", "submit.sh"],
                cwd=prepack_dir,
                capture_output=True,
                text=True,
                check=True,
            )
            job_id = result.stdout.strip().split()[-1]
            print(f"  OK     {label}  — job {job_id}")
            submitted.append((label, job_id))

        except subprocess.CalledProcessError as e:
            msg = e.stderr.strip() or e.stdout.strip()
            print(f"  ERROR  {label}  — {msg}")
            failed.append(label)

        except FileNotFoundError:
            print("Error: sbatch not found. Run this script on the cluster.")
            sys.exit(1)

    # ── Summary ───────────────────────────────────────────────────────────────
    print()
    print(f"Submitted : {len(submitted)}")
    print(f"Skipped   : {len(skipped)}  (already done)")
    print(f"Failed    : {len(failed)}")

    if submitted:
        job_ids = ",".join(j for _, j in submitted)
        print()
        print("Useful commands:")
        print(f"  Watch queue   : squeue -u $USER")
        print(f"  Cancel all    : scancel {job_ids}")
        print(f"  Check logs    : tail -f <GH>/<LIG>/1prepack/logs/prepack_<jobid>.out")


if __name__ == "__main__":
    main()
