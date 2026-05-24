#!/usr/bin/env python3
"""
glycandock_scores.py
--------------------
Processes GlycanDock results across multiple ligand subdirectories.
Run from a GH folder (e.g., cd 39) to score all oligosaccharides at once.

For each [LIG]/2dock/ subfolder with a results/ folder:
  1. Loads all score_*.sc files (any nstruct, any number of files)
  2. Filters: n_tor_cycles >= TOR_CYCLE_MIN  (torsion sampling must have run)
              ring_Lrmsd < RING_LRMSD_CUTOFF  (glycan must stay near binding site)
  3. Top 20% by total_score
  4. Top N by interaction_energy
  5. Writes CSV summary + full scores CSV + TXT report
  6. Copies top PDB files (handles both .pdb and .pdb.gz) to [LIG]/2dock/Top_PDBs/

Usage: cd /path/to/GlycanDock/39 && python3 glycandock_scores.py
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


def find_pdb_file(results_path, description):
    """
    Find the .pdb.gz or .pdb file for a given description.
    Returns (path, is_gz) or (None, None) if not found.
    """
    for ext, is_gz in [('.pdb.gz', True), ('.pdb', False)]:
        name = description if description.endswith(ext) else description + ext
        path = os.path.join(results_path, name)
        if os.path.exists(path):
            return path, is_gz
    return None, None


def copy_pdb(results_path, dst_dir, description):
    """Decompress .pdb.gz (or copy .pdb) into dst_dir as description.pdb."""
    src_path, is_gz = find_pdb_file(results_path, description)
    if src_path is None:
        return None

    # Destination is always plain .pdb
    base = description
    if base.endswith('.pdb.gz'):
        base = base[:-7]
    elif base.endswith('.pdb'):
        base = base[:-4]
    dst_path = os.path.join(dst_dir, base + '.pdb')

    if is_gz:
        with gzip.open(src_path, 'rb') as f_in, open(dst_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    else:
        shutil.copy2(src_path, dst_path)

    return dst_path


def process_subdir(dock_path, lig_name):
    """
    Score and filter GlycanDock results for one ligand's 2dock/ directory.
    dock_path = [GH]/[LIG]/2dock/
    """
    results_path = os.path.join(dock_path, 'results')
    top_pdbs_path = os.path.join(dock_path, 'Top_PDBs')

    if not os.path.isdir(results_path):
        print(f"  ⚠️  No results/ folder in {lig_name}/2dock/ — skipping.")
        return

    if os.path.exists(top_pdbs_path):
        print(f"  ⏭️  Top_PDBs/ already exists — skipping.")
        return

    # ── Load all score files ──────────────────────────────────────────────────
    rows, header = load_scores_from_dir(results_path)
    if not rows:
        print(f"  ❌ No score data found in {lig_name}/2dock/results/")
        return

    n_total = len(rows)
    print(f"  📊 Loaded {n_total} structures")

    # ── Stage 1: constraint filter ────────────────────────────────────────────
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
        f.write(f"GlycanDock Top {len(top)} Structures — {lig_name}\n")
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

    print(f"  📁 Copied {len(copied)} PDB(s) to {lig_name}/2dock/Top_PDBs/")
    print(f"  📝 Report:  {txt_out}")
    print(f"  📝 CSV:     {csv_summary}")


def main():
    base_dir = os.getcwd()

    # Find all [LIG]/2dock/ subdirectories one level down
    dock_dirs = []
    for lig in sorted(os.listdir(base_dir)):
        if lig.startswith('.'):
            continue
        lig_path = os.path.join(base_dir, lig)
        if not os.path.isdir(lig_path):
            continue
        dock_path = os.path.join(lig_path, '2dock')
        if os.path.isdir(dock_path):
            dock_dirs.append((lig, dock_path))

    if not dock_dirs:
        print("No [LIG]/2dock/ directories found.")
        print("Run from a GH folder, e.g.:  cd 39 && python3 glycandock_scores.py")
        return

    print(f"Found {len(dock_dirs)} ligand(s) with 2dock/ in {base_dir}\n")
    for lig_name, dock_path in dock_dirs:
        print(f"📂 {lig_name}/2dock")
        process_subdir(dock_path, lig_name)
        print()


if __name__ == '__main__':
    main()
