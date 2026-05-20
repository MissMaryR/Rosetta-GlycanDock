#!/usr/bin/env python3
"""
glycandock_scores.py
--------------------
Processes GlycanDock results across multiple subdirectories.
For each subdir with a results/ folder:
  1. Loads all score_*.sc files (any nstruct, any number of files)
  2. Filters: n_tor_cycles > 0  (torsion sampling must have run)
              ring_Lrmsd < RING_LRMSD_CUTOFF  (glycan must stay near binding site)
  3. Sorts by interaction_energy (ascending — most negative = best binding)
  4. Selects top N structures
  5. Writes CSV summary + full scores CSV + TXT report
  6. Copies top PDB files (handles both .pdb and .pdb.gz) to Top_PDBs/

Usage: cd /path/to/GlycanDock/76 && python3 script/glycandock_scores.py
"""

import os
import csv
import gzip
import math
import shutil  # used for copyfileobj

# ── Settings ──────────────────────────────────────────────────────────────────
TOP_N             = 10   # number of top structures to select
TOP_PCT           = 0.20 # fraction kept by total_score before final interaction_energy sort
TOR_CYCLE_MIN     = 1    # n_tor_cycles must be >= this (confirms torsion sampling ran)
RING_LRMSD_CUTOFF = 10.0 # Å — exclude structures where glycan drifted out of binding site

# GlycanDock-specific fields to include in the summary CSV (in display order)
SUMMARY_FIELDS = [
    'description',
    'interaction_energy',   # primary ranking metric (most negative = best)
    'ring_Lrmsd',           # glycan displacement from starting pose
    'ring_Srmsd',           # glycan RMSD from starting pose (all atoms)
    'n_tor_cycles',         # torsion sampling cycles (must be > 0)
    'n_tor_moves_accepted', # torsion moves accepted
    'n_rb_cycles',          # rigid-body sampling cycles
    'n_rb_moves_accepted',  # rigid-body moves accepted
    'n_intf_residues',      # interface residues in model
    'n_intf_res_contacts',  # interface residue-residue contacts
    'Fnat',                 # fraction native contacts (vs input pose)
    'Fnat_intf_residues',   # fraction native interface residues (vs input pose)
    'sugar_bb',             # glycosidic torsion energy
    'total_score',          # overall Rosetta energy (not used for ranking)
]
# ──────────────────────────────────────────────────────────────────────────────


def safe_float(val):
    try:
        return float(val)
    except (ValueError, TypeError):
        return float('nan')


def load_scores_from_dir(results_path):
    """
    Load all GlycanDock score entries from a results/ directory.
    Handles any number of score_*.sc files and any nstruct per file.

    Rosetta GlycanDock score files have this format:
        SEQUENCE:
        SCORE: total_score  Fnat  ...  description   <- header (non-numeric after SCORE:)
        SCORE:  -2263.891   0.721 ...  input_0001    <- data rows

    We skip SEQUENCE: lines, strip the leading 'SCORE:' token, treat the first
    SCORE: line (non-numeric first value) as the column header, and subsequent
    SCORE: lines (numeric first value) as data rows.
    """
    rows = []
    header = []

    for filename in sorted(os.listdir(results_path)):
        if not (filename.startswith('score') and filename.endswith('.sc')):
            continue
        filepath = os.path.join(results_path, filename)
        local_header = []

        with open(filepath) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('SEQUENCE'):
                    continue
                tokens = line.split()
                # Only process SCORE: lines; strip the leading 'SCORE:' token
                if not tokens or tokens[0] != 'SCORE:':
                    continue
                tokens = tokens[1:]
                if not tokens:
                    continue

                if not local_header:
                    # First SCORE: line has column names (first token non-numeric)
                    local_header = tokens
                    if not header:
                        header = list(tokens)
                    continue

                # Skip any repeated header lines
                try:
                    float(tokens[0])
                except ValueError:
                    continue

                if len(tokens) < len(local_header):
                    continue

                row = {'_results_path': results_path}
                for i, col in enumerate(local_header):
                    if col == 'description':
                        row[col] = tokens[-1].strip()
                    else:
                        row[col] = safe_float(tokens[i])
                rows.append(row)

    return rows, header


def find_pdb_gz(results_path, description):
    """
    Find the .pdb.gz file for a given description.
    Returns path or None if not found.
    """
    name = description if description.endswith('.pdb.gz') else description + '.pdb.gz'
    path = os.path.join(results_path, name)
    return path if os.path.exists(path) else None


def copy_pdb(results_path, dst_dir, description):
    """Decompress .pdb.gz into dst_dir as description.pdb."""
    src_path = find_pdb_gz(results_path, description)
    if src_path is None:
        return None

    dst_name = (description if description.endswith('.pdb') else description + '.pdb')
    dst_path = os.path.join(dst_dir, dst_name)

    with gzip.open(src_path, 'rb') as f_in, open(dst_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

    return dst_path


def process_subdir(subdir_path, subdir_name):
    results_path = os.path.join(subdir_path, 'results')
    top_pdbs_path = os.path.join(subdir_path, 'Top_PDBs')

    if not os.path.isdir(results_path):
        print(f"  ⚠️  No results/ folder — skipping.")
        return

    if os.path.exists(top_pdbs_path):
        print(f"  ⏭️  Top_PDBs/ already exists — skipping.")
        return

    # ── Load all score files ──────────────────────────────────────────────────
    rows, header = load_scores_from_dir(results_path)
    if not rows:
        print(f"  ❌ No score data found in results/")
        return

    n_total = len(rows)
    print(f"  📊 Loaded {n_total} structures")

    # ── Stage 1: constraint filter ────────────────────────────────────────────
    # Two conditions (equivalent to rosetta_scores9.py's all_cst < 1.0):
    #   1. n_tor_cycles >= TOR_CYCLE_MIN  — torsion sampling ran (invalid if 0)
    #   2. ring_Lrmsd < RING_LRMSD_CUTOFF — glycan stayed near the binding site
    filtered = [
        r for r in rows
        if r.get('n_tor_cycles', 0) >= TOR_CYCLE_MIN
        and r.get('ring_Lrmsd', float('inf')) < RING_LRMSD_CUTOFF
    ]
    n_filtered = len(filtered)
    print(f"  ✅ Passed constraint filter (n_tor_cycles >= {TOR_CYCLE_MIN} and "
          f"ring_Lrmsd < {RING_LRMSD_CUTOFF} Å): {n_filtered} / {n_total}")
    if not filtered:
        print(f"  ❌ No structures passed constraint filter.")
        return

    # ── Stage 2: top 20% by total_score ──────────────────────────────────────
    filtered.sort(key=lambda r: r.get('total_score', float('inf')))
    n_top20 = math.ceil(TOP_PCT * len(filtered))
    top_20 = filtered[:n_top20]
    print(f"  ✅ Top 20% by total_score: {len(top_20)}")

    # ── Stage 3: top N by interaction_energy ─────────────────────────────────
    top_20.sort(key=lambda r: r.get('interaction_energy', float('inf')))
    top = top_20[:TOP_N]
    print(f"  🏆 Top {len(top)} by interaction_energy: "
          f"{top[0].get('interaction_energy', float('nan')):.3f} to "
          f"{top[-1].get('interaction_energy', float('nan')):.3f} REU")

    # ── Write outputs ─────────────────────────────────────────────────────────
    os.makedirs(top_pdbs_path, exist_ok=True)

    available_summary = [f for f in SUMMARY_FIELDS if f in header or f == 'description']

    def clean_row(row, fields):
        return {k: row.get(k, 'NA') for k in fields}

    # Summary CSV (key GlycanDock fields only)
    csv_summary = os.path.join(top_pdbs_path, 'top_glycandock_summary.csv')
    with open(csv_summary, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=available_summary)
        writer.writeheader()
        for row in top:
            writer.writerow(clean_row(row, available_summary))

    # Full scores CSV (all Rosetta columns)
    csv_full = os.path.join(top_pdbs_path, 'top_glycandock_fullscores.csv')
    full_fields = [c for c in header if c != 'SCORE:']
    with open(csv_full, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=full_fields)
        writer.writeheader()
        for row in top:
            writer.writerow(clean_row(row, full_fields))

    # Human-readable TXT report
    txt_out = os.path.join(top_pdbs_path, 'top_glycandock_report.txt')
    col_w = 22
    with open(txt_out, 'w') as f:
        f.write(f"GlycanDock Top {len(top)} Structures — {subdir_name}\n")
        f.write("=" * 80 + "\n")
        f.write(f"Filters:  n_tor_cycles >= {TOR_CYCLE_MIN} and ring_Lrmsd < {RING_LRMSD_CUTOFF} Å"
                f"  →  top {int(TOP_PCT*100)}% by total_score  →  top {TOP_N} by interaction_energy\n\n")

        header_line = '  '.join(f"{k:<{col_w}}" for k in available_summary)
        f.write(header_line + "\n")
        f.write("-" * len(header_line) + "\n")

        for row in top:
            parts = []
            for k in available_summary:
                v = row.get(k, 'NA')
                if isinstance(v, float):
                    parts.append(f"{v:<{col_w}.3f}")
                else:
                    parts.append(f"{str(v):<{col_w}}")
            f.write('  '.join(parts) + "\n")

        f.write("\n" + "=" * 80 + "\n")
        f.write("Summary Statistics:\n")
        f.write(f"  Total structures loaded:            {n_total}\n")
        f.write(f"  Passed constraint filter:           {n_filtered}\n"
                f"    (n_tor_cycles >= {TOR_CYCLE_MIN} and ring_Lrmsd < {RING_LRMSD_CUTOFF} Å)\n")
        f.write(f"  Top 20% by total_score:             {len(top_20)}\n")
        f.write(f"  Selected top {TOP_N} by interaction_energy: {len(top)}\n")
        f.write(f"  Best interaction_energy:            {top[0].get('interaction_energy', float('nan')):.3f} REU\n")
        f.write(f"  Worst of top {TOP_N}:                    {top[-1].get('interaction_energy', float('nan')):.3f} REU\n")

    # Copy + decompress PDB files
    copied, missing = [], []
    for row in top:
        desc = row.get('description', '')
        dst = copy_pdb(results_path, top_pdbs_path, desc)
        if dst:
            copied.append(dst)
        else:
            missing.append(desc)

    if missing:
        for m in missing:
            print(f"  ⚠️  PDB not found: {m}")

    print(f"  📁 Copied {len(copied)} PDB(s) to Top_PDBs/")
    print(f"  📝 Report:  {txt_out}")
    print(f"  📝 CSV:     {csv_summary}")


def main():
    base_dir = os.getcwd()
    subdirs = sorted(
        d for d in os.listdir(base_dir)
        if os.path.isdir(os.path.join(base_dir, d)) and not d.startswith('.')
    )

    if not subdirs:
        print("No subdirectories found. Run from the GlycanDock target directory.")
        return

    print(f"Found {len(subdirs)} subdirectory/ies in {base_dir}\n")
    for name in subdirs:
        path = os.path.join(base_dir, name)
        print(f"📂 {name}")
        process_subdir(path, name)
        print()


if __name__ == '__main__':
    main()
