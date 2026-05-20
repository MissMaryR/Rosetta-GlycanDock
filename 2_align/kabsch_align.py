"""
kabsch_align.py — PyMOL-free alternative to align_saccharide.pml

Replicates PyMOL's pair_fit using the Kabsch algorithm (SVD-based rigid-body
superposition). Reads docked.pdb and [LIG].pdb from the current directory,
writes input_aligned.pdb ready for GlycanDock prepack.

Requirements: Python 3, NumPy

Usage:
    python kabsch_align.py CL3     # or CR3, H2B, H3B, XY3

The ligand name selects the atom mapping and HETNAM/LINK records.
Both docked.pdb and [LIG].pdb must be in the same directory.

Atom mappings
─────────────
pair_fit aligns the non-reducing-end ring of the Rosetta saccharide (mobile)
onto the deepest ring in the docked pose (target) using 6 atom pairs.
Fitting only ONE ring ensures a pure rigid-body transform that preserves 4C1
chair geometry. GlycanDock samples the other rings' glycosidic torsions in
Stage 1+2.

CL3/CR3/H2B/H3B (glucose ring):
  Rosetta C1 (anomeric) ↔ docked C4   Rosetta C4 ↔ docked C1
  Rosetta C2            ↔ docked C2   Rosetta C5 ↔ docked C5
  Rosetta C3            ↔ docked C3   Rosetta O5 ↔ docked O1

XY3 (xylose ring — scrambled atom order in docked pose):
  Rosetta C1 ↔ docked C1   Rosetta C4 ↔ docked C4
  Rosetta C2 ↔ docked C5   Rosetta C5 ↔ docked C3
  Rosetta C3 ↔ docked C2   Rosetta O5 ↔ docked O1
"""

import sys
import numpy as np
from pathlib import Path

# ── Per-ligand config ─────────────────────────────────────────────────────────
LIGANDS = {
    "CL3": {
        "mobile_resi":  3,
        "mobile_atoms": ["C1", "C2", "C3", "C4", "C5", "O5"],
        "target_resn":  "CL3",
        "target_atoms": ["C4", "C2", "C3", "C1", "C5", "O1"],
        "hetnam_lines": [
            "HETNAM     Glc X   1  ->4)-beta-D-Glcp\n",
            "HETNAM     Glc X   2  ->4)-beta-D-Glcp\n",
            "HETNAM     Glc X   3  ->4)-beta-D-Glcp\n",
        ],
        "link_lines": [
            "LINK         O4  Glc X   1                 C1  Glc X   2                  1.50  \n",
            "LINK         O4  Glc X   2                 C1  Glc X   3                  1.50  \n",
        ],
    },
    "CR3": {
        "mobile_resi":  3,
        "mobile_atoms": ["C1", "C2", "C3", "C4", "C5", "O5"],
        "target_resn":  "CR3",
        "target_atoms": ["C4", "C2", "C3", "C1", "C5", "O1"],
        "hetnam_lines": [
            "HETNAM     Glc X   1  ->3)-beta-D-Glcp\n",
            "HETNAM     Glc X   2  ->3)-beta-D-Glcp\n",
            "HETNAM     Glc X   3  ->4)-beta-D-Glcp\n",
        ],
        "link_lines": [
            "LINK         O3  Glc X   1                 C1  Glc X   2                  1.50  \n",
            "LINK         O3  Glc X   2                 C1  Glc X   3                  1.50  \n",
        ],
    },
    "H2B": {
        "mobile_resi":  2,
        "mobile_atoms": ["C1", "C2", "C3", "C4", "C5", "O5"],
        "target_resn":  "H2B",
        "target_atoms": ["C4", "C2", "C3", "C1", "C5", "O1"],
        "hetnam_lines": [
            "HETNAM     Glc X   1  ->4)-beta-D-Glcp\n",
            "HETNAM     Glc X   2  ->4)-beta-D-Glcp\n",
            "HETNAM     Xyl X   3  ->4)-alpha-D-Xylp\n",
        ],
        "link_lines": [
            "LINK         O4  Glc X   1                 C1  Glc X   2                  1.50  \n",
            "LINK         O6  Glc X   1                 C1  Xyl X   3                  1.50  \n",
        ],
    },
    "H3B": {
        "mobile_resi":  3,
        "mobile_atoms": ["C1", "C2", "C3", "C4", "C5", "O5"],
        "target_resn":  "H3B",
        "target_atoms": ["C4", "C2", "C3", "C1", "C5", "O1"],
        "hetnam_lines": [
            "HETNAM     Glc X   1  ->4)-beta-D-Glcp\n",
            "HETNAM     Glc X   2  ->4)-beta-D-Glcp\n",
            "HETNAM     Glc X   3  ->4)-beta-D-Glcp\n",
            "HETNAM     Xyl X   4  ->4)-alpha-D-Xylp\n",
        ],
        "link_lines": [
            "LINK         O4  Glc X   1                 C1  Glc X   2                  1.50  \n",
            "LINK         O4  Glc X   2                 C1  Glc X   3                  1.50  \n",
            "LINK         O6  Glc X   3                 C1  Xyl X   4                  1.50  \n",
        ],
    },
    "XY3": {
        "mobile_resi":  3,
        "mobile_atoms": ["C1", "C2", "C3", "C4", "C5", "O5"],
        "target_resn":  "XY3",
        "target_atoms": ["C1", "C5", "C2", "C4", "C3", "O1"],
        "hetnam_lines": [
            "HETNAM     Xyl X   1  ->4)-beta-D-Xylp\n",
            "HETNAM     Xyl X   2  ->4)-beta-D-Xylp\n",
            "HETNAM     Xyl X   3  ->4)-beta-D-Xylp\n",
        ],
        "link_lines": [
            "LINK         O4  Xyl X   1                 C1  Xyl X   2                  1.50  \n",
            "LINK         O4  Xyl X   2                 C1  Xyl X   3                  1.50  \n",
        ],
    },
}

# ── PDB parsing ───────────────────────────────────────────────────────────────

def parse_atoms(path):
    atoms = []
    with open(path) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            atoms.append({
                "record": line[:6].strip(),
                "name":   line[12:16].strip(),
                "resn":   line[17:20].strip(),
                "resi":   int(line[22:26]),
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
                "raw":    line,
            })
    return atoms


def get_coords(atoms, names, resi=None, resn=None):
    sel = atoms
    if resi is not None: sel = [a for a in sel if a["resi"] == resi]
    if resn is not None: sel = [a for a in sel if a["resn"] == resn]
    lookup = {a["name"]: np.array([a["x"], a["y"], a["z"]]) for a in sel}
    missing = [n for n in names if n not in lookup]
    if missing:
        raise ValueError(
            f"Atoms not found: {missing}  "
            f"(resi={resi}, resn={resn})  "
            f"Available: {sorted(lookup.keys())}"
        )
    return np.array([lookup[n] for n in names])


def kabsch(P, Q):
    """Optimal rotation of P onto Q. Returns R, t such that Q ≈ R @ P + t."""
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, Qc - R @ Pc


def transform(atoms, R, t):
    out = []
    for a in atoms:
        p = np.array([a["x"], a["y"], a["z"]])
        q = R @ p + t
        n = dict(a)
        n["x"], n["y"], n["z"] = float(q[0]), float(q[1]), float(q[2])
        out.append(n)
    return out


def fmt_line(a, serial):
    raw  = a["raw"]
    name = a["name"]
    if   len(name) >= 4: nf = name[:4]
    elif len(name) == 1: nf = f" {name}  "
    elif len(name) == 2: nf = f" {name} "
    else:                nf = f" {name}"
    rec = "HETATM" if a["record"] == "HETATM" else "ATOM  "
    line = (
        f"{rec}{serial:5d} {nf}{raw[16]}{raw[17:20]} {raw[21]}{raw[22:26]}{raw[26]}   "
        f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}"
        f"{raw[54:60] if len(raw) > 54 else '  1.00'}"
        f"{raw[60:66] if len(raw) > 60 else '  0.00'}"
        f"{raw[66:]   if len(raw) > 66 else ''}"
    )
    return line.rstrip("\n") + "\n"

# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    if len(sys.argv) < 2 or sys.argv[1] not in LIGANDS:
        print(f"Usage: python kabsch_align.py <LIG>")
        print(f"  LIG must be one of: {', '.join(LIGANDS)}")
        sys.exit(1)

    lig = sys.argv[1]
    cfg = LIGANDS[lig]
    here = Path(__file__).resolve().parent

    docked_path = here / "docked.pdb"
    sac_path    = here / f"{lig}.pdb"
    out_path    = here / "input_aligned.pdb"

    for p in (docked_path, sac_path):
        if not p.exists():
            print(f"Error: {p} not found")
            sys.exit(1)

    docked_atoms = parse_atoms(docked_path)
    sac_atoms    = parse_atoms(sac_path)

    target_coords = get_coords(docked_atoms, cfg["target_atoms"], resn=cfg["target_resn"])
    mobile_coords = get_coords(sac_atoms,    cfg["mobile_atoms"], resi=cfg["mobile_resi"])

    R, t = kabsch(mobile_coords, target_coords)

    # RMSD on the 6 fitted atoms
    fitted = (R @ mobile_coords.T).T + t
    rmsd   = float(np.sqrt(np.mean(np.sum((fitted - target_coords) ** 2, axis=1))))
    print(f"Pair-fit RMSD ({lig}): {rmsd:.4f} Å")
    if rmsd > 0.5:
        print("  WARNING: RMSD > 0.5 Å — verify atom mapping")

    # Apply transform to full saccharide
    sac_aligned   = transform(sac_atoms, R, t)
    protein_atoms = [a for a in docked_atoms if a["record"] == "ATOM"]
    glycan_atoms  = [a for a in sac_aligned  if a["record"] == "HETATM"]
    if not glycan_atoms:
        glycan_atoms = sac_aligned  # fallback if saccharide uses ATOM records

    lines = list(cfg["hetnam_lines"]) + list(cfg["link_lines"])
    serial = 1
    for a in protein_atoms:
        lines.append(fmt_line(a, serial)); serial += 1
    for a in glycan_atoms:
        a2 = dict(a); a2["record"] = "HETATM"
        lines.append(fmt_line(a2, serial)); serial += 1
    lines += ["TER\n", "END\n"]

    with open(out_path, "w") as f:
        f.writelines(lines)

    print(f"Written: {out_path}")
    print(f"  {len(protein_atoms)} ATOM + {len(glycan_atoms)} HETATM records")
    print("Copy input_aligned.pdb to 3_prepack/ and run sbatch submit.sh")


if __name__ == "__main__":
    main()
