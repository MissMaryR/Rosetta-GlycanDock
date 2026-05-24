"""
GlycanDock Docking Submission Script

Finds all [LIG]/2dock/ directories in the same folder as this script,
skips any that already have results, and submits the rest via the
submit.sh present in each 2dock/ folder.

Usage:
    python submit_dock.py          # submits all ligands for this GH
"""

import sys
import subprocess
from pathlib import Path


def main():
    # Default to the directory containing this script so it works
    # correctly regardless of where it is called from
    root = Path(__file__).resolve().parent

    dock_dirs = sorted(root.glob("*/2dock"))

    if not dock_dirs:
        print(f"No [LIG]/2dock/ directories found under {root}")
        sys.exit(1)

    gh = root.name
    print(f"GH: {gh}")
    print(f"Dock folders: {len(dock_dirs)} found")
    print()

    submitted = []
    skipped   = []
    failed    = []

    for dock_dir in dock_dirs:
        lig = dock_dir.parent.name
        label = f"{gh}/{lig}"

        submit_sh  = dock_dir / "submit.sh"
        results    = dock_dir / "results"

        # Skip if results directory already has any output files
        if any(results.glob("*.pdb.gz")) or any(results.glob("*.pdb")):
            print(f"  SKIP   {label}  (results already exist)")
            skipped.append(label)
            continue

        if not submit_sh.exists():
            print(f"  ERROR  {label}  (no submit.sh found)")
            failed.append(label)
            continue

        if not (dock_dir / "packed.pdb").exists():
            print(f"  ERROR  {label}  (no packed.pdb found)")
            failed.append(label)
            continue

        try:
            result = subprocess.run(
                ["sbatch", "submit.sh"],
                cwd=dock_dir,
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
        print(f"  Watch queue  : squeue -u $USER")
        print(f"  Cancel all   : scancel {job_ids}")
        print(f"  Check a log  : tail -f {gh}/<LIG>/2dock/logs/glycandock_<A>_<a>.out")


if __name__ == "__main__":
    main()
