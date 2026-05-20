# PyMOL script: align Rosetta-generated saccharide onto a reference docked pose
# Step 2 in the GlycanDock pipeline.
#
# Before running:
#   1. Copy out_saccharide.pdb (from 1_gen_saccharide/) to this directory
#   2. Copy your docked.pdb (relaxed protein + reference ligand) to this directory
#   3. Edit the pair_fit atom names below to match your reference ligand
#
# Run from PyMOL: File > Run Script > align_saccharide.pml
# OR from terminal: pymol -c align_saccharide.pml
#
# Output: input_aligned.pdb — copy this to 3_prepack/ as input_aligned.pdb

load docked.pdb,         docked      # relaxed protein (chain A) + reference ligand
load out_saccharide.pdb, saccharide  # Rosetta-native saccharide from 1_gen_saccharide/

# ── Align Glc resi 3 (non-reducing end) onto CL3 ring 1 — 6 atoms ONLY ───────
# CRITICAL: fit only ONE ring (6 atom pairs). Fitting all 3 rings simultaneously
# (18 pairs) forces a global compromise that destroys individual ring 4C1 chair
# geometry — Rosetta then misidentifies the rings as furanose/wrong linkage type
# and crashes (->2)-alpha-D-Glcf instead of ->4)-beta-D-Glcp).
#
# With 6 pairs, pair_fit is a true rigid-body transform: all atoms including
# rings 2 and 1 move with the same rotation, preserving Rosetta's native ring
# geometry. GlycanDock Stage 1+2 will adjust the glycosidic torsions of rings
# 2 and 1 during sampling.
#
# Atom mapping (non-reducing end Glc C1↔CL3 C4 because CL3 numbers in the
# reducing→non-reducing direction, Glc numbers anomeric-first):
#   Glc C1 = anomeric C  ↔  CL3 C4 (glycosidic position in ring 1)
#   Glc C4 = glycosidic C ↔  CL3 C1 (anomeric position in ring 1)
#   Glc O5 = ring O       ↔  CL3 O1 (ring O in CL3 numbering)

pair_fit \
  saccharide and resi 3 and name C1,  docked and resn CL3 and name C4,  \
  saccharide and resi 3 and name C2,  docked and resn CL3 and name C2,  \
  saccharide and resi 3 and name C3,  docked and resn CL3 and name C3,  \
  saccharide and resi 3 and name C4,  docked and resn CL3 and name C1,  \
  saccharide and resi 3 and name C5,  docked and resn CL3 and name C5,  \
  saccharide and resi 3 and name O5,  docked and resn CL3 and name O1

# ── Extract protein and write final PDB with correct section order ────────────
# Rosetta requires the glycoligand coordinates at the BOTTOM of the PDB file
# (after all protein ATOM records) for correct FoldTree construction.
# PyMOL's save command puts HETATM before ATOM in combined objects, so we
# handle the merge and ordering explicitly in the Python block below.
create protein, docked and chain A

python
import re

out_path  = "input_aligned.pdb"  # output written to current directory
prot_tmp  = "/tmp/glycandock_prot.pdb"
glc_tmp   = "/tmp/glycandock_glc.pdb"

cmd.save(prot_tmp, "protein")
cmd.save(glc_tmp,  "saccharide")

hetnam_lines = [
    "HETNAM     Glc X   1  ->4)-beta-D-Glcp\n",
    "HETNAM     Glc X   2  ->4)-beta-D-Glcp\n",
    "HETNAM     Glc X   3  ->4)-beta-D-Glcp\n",
]
# LINK records declare the two β-1→4 glycosidic bonds (O4 of residue n → C1 of residue n+1).
# These match the LINK records in out_saccharide.pdb exactly. Required for -maintain_links.
# DO NOT use -auto_detect_glycan_connections: distance scan finds a spurious 3rd linkage
# (POLYMER_LOWER=0 → "Found 2 glycan trees" → crash). Use -maintain_links in flags instead.
link_lines = [
    "LINK         O4  Glc X   1                 C1  Glc X   2                  1.50  \n",
    "LINK         O4  Glc X   2                 C1  Glc X   3                  1.50  \n",
]

with open(prot_tmp) as f:
    atom_lines = [l for l in f if l.startswith("ATOM")]
with open(glc_tmp) as f:
    hetatm_lines = [l for l in f if l.startswith("HETATM")]

# Write: HETNAM -> protein ATOM -> glycan HETATM -> TER -> END
# HETNAM records set the correct residue type (->4)-beta-D-Glcp), preventing
# the Glcf/alpha/2-linked misassignment crash.
# NO LINK records: when both LINK records and -auto_detect_glycan_connections
# are active, Rosetta counts the root-residue glycosidic bond twice (once from
# LINK, once from distance scan) → 3 linkages instead of 2 → broken FoldTree
# → POLYMER_LOWER=0 crash. Let auto_detect be the sole connectivity source.
# NO CONECT records: partial CONECT (inter-residue only) similarly double-counts.
result = list(hetnam_lines) + list(link_lines)
serial = 1
for line in atom_lines:
    result.append(line[:6] + f"{serial:5d}" + line[11:])
    serial += 1
for line in hetatm_lines:
    result.append(line[:6] + f"{serial:5d}" + line[11:])
    serial += 1
# No CONECT records needed: LINK records + -maintain_links handle all connectivity.
result.append("TER\n")
result.append("END\n")

with open(out_path, "w") as f:
    f.writelines(result)
print(f"Written {len(atom_lines)} ATOM + {len(hetatm_lines)} HETATM records")
print(f"  + 3 HETNAM + 2 LINK records in header")
print("Glycan is at the bottom — correct for Rosetta FoldTree")
print("Use -maintain_links (not -auto_detect_glycan_connections) in flags")
print("Done — input_aligned.pdb ready. Copy to cluster as input.pdb")
python end
