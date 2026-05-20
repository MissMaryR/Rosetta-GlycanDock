# GlycanDock Pipeline for Oligosaccharide Docking

A complete working pipeline for docking oligosaccharides (DP ≥ 2) into protein binding sites using **GlycanDock** in Rosetta.

Based on: Nance et al. 2021, *J. Chem. Theory Comput.* — [GlycanDock paper](https://doi.org/10.1021/acs.jctc.1c00650).

---

## Overview

GlycanDock is a standalone Rosetta binary for docking oligosaccharides. It samples glycosidic torsions and rigid-body orientation simultaneously using Monte Carlo + minimization with soft-to-hard van der Waals ramping. **Not for single-residue sugars** — use GALigandDock for those.

```
1_gen_saccharide/   →   2_align/   →   3_prepack/   →   4_dock/   →   scripts/
  cluster: ~10s         local          cluster: ~2min    cluster:        local or cluster
  build saccharide       align to       prepack once      array job       score & rank
  from sequence          active site    (per GH/ligand)   (≥200 structs)
```

### Recommended folder structure for multi-target runs

When docking multiple oligosaccharides against multiple GH enzymes, organise each (GH, ligand) pair as:

```
[GH]/
  [LIG]/
    0align/     ← docked.pdb, [LIG].pdb, align_saccharide.pml / kabsch_align.py, input_aligned.pdb
    1prepack/   ← flags, submit.sh, input_aligned.pdb → results/packed.pdb
    2dock/      ← flags, submit.sh, packed.pdb → results/
  submit_dock.py   ← submits all 2dock/ jobs for this GH (one copy per GH folder)
scripts/
  submit_all_prepack.py   ← submit all prepack jobs across all GH/LIG pairs at once
  submit_dock.py          ← template; copy one into each [GH]/ folder
  glycandock_scores.py    ← score, filter, and rank docking results
  Scoring.sh              ← SLURM wrapper to run glycandock_scores.py on the cluster
```

---

## Prerequisites

- **Rosetta** with carbohydrate support: `GlycanDock.static.linuxgccrelease`, `pose_from_saccharide_sequence.static.linuxgccrelease`
- **PyMOL** (any version with `pair_fit`) *or* **Python 3 + NumPy** (for `kabsch_align.py`)
- **SLURM** cluster (or equivalent job scheduler)
- **Python 3** for submission scripts and `glycandock_scores.py` (no non-stdlib dependencies beyond NumPy for alignment)
- **Protein input**: Rosetta-relaxed structure, chain A (e.g., AlphaFold3 model run through Rosetta FastRelax)

---

## Step 1 — Generate Rosetta-Native Saccharide (`1_gen_saccharide/`)

**Run on cluster. Takes ~10 seconds.**

The `saccharides/` folder contains ready-to-use scripts and PDB outputs for each ligand. To generate a saccharide, submit the appropriate script:

```bash
cd 1_gen_saccharide/saccharides/
bash CL3.sh     # generates CL3.pdb  (cellotriose, β-1→4 Glc ×3)
bash CR3.sh     # generates CR3.pdb  (laminaritrioside, β-1→3 Glc ×3)
bash H2B.sh     # generates H2B.pdb  (Glc-Glc + α-Xyl-1,6 on reducing Glc)
bash H3B.sh     # generates H3B.pdb  (Glc-Glc-Glc + α-Xyl-1,6 on non-reducing Glc)
bash XY3.sh     # generates XY3.pdb  (xylotriose, β-1→4 Xyl ×3)
```

Each script calls `pose_from_saccharide_sequence` with the correct sequence string and renames the output to `[LIG].pdb`. Copy `[LIG].pdb` into your `[GH]/[LIG]/0align/` directory for Step 2.

To define a new saccharide, edit the `-carbohydrates:saccharide_sequence` string in `gen_saccharide.sh` (see [Saccharide Sequence Notation](#saccharide-sequence-notation) below).

> **Why not use a PDB ligand directly?** Manually converting a ligand PDB to Rosetta residue names always produces wrong ring geometry → Rosetta assigns the wrong residue type (e.g., `->3)-beta-D-Glcp` instead of `->4)-beta-D-Glcp`) → torsion sampling never activates. Always build the saccharide natively with `pose_from_saccharide_sequence`.

---

## Step 2 — Align Saccharide into Binding Site (`2_align/`)

**Run locally. Two options: PyMOL script or Python (no PyMOL required).**

You need two inputs in the same working directory:
- `[LIG].pdb` — Rosetta-native saccharide from Step 1
- `docked.pdb` — Rosetta-relaxed protein (chain A) with a GALigandDock reference pose of your ligand (chain X, resi 1)

### Option A — PyMOL script

Edit `align_saccharide.pml` to match your reference ligand's atom naming, then run:

```bash
pymol -c align_saccharide.pml
# OR in PyMOL GUI: File > Run Script > align_saccharide.pml
```

### Option B — Python (NumPy only, no PyMOL)

```bash
python kabsch_align.py CL3    # or CR3, H2B, H3B, XY3
```

Both options align the **non-reducing end ring** of the Rosetta saccharide onto the deepest ring in the active site using 6 atom pairs, then write `input_aligned.pdb`.

> **Critical: fit only ONE ring (6 atom pairs).** Fitting all rings simultaneously distorts 4C1 chair geometry → Rosetta misidentifies ring types and crashes. With 6 pairs, the superposition is a pure rigid-body transform — all atoms move together preserving native ring geometry. GlycanDock samples glycosidic torsions of the other rings during docking.

Copy `input_aligned.pdb` to `3_prepack/`.

---

## Step 3 — Prepack (`3_prepack/`)

**Run on cluster once per (GH, ligand) pair. Takes ~2 minutes.**

The prepack step repacks protein side chains and glycan hydroxyls at the interface. This erases alignment bias and gives each docking trajectory a clean starting point.

```bash
cp ../0align/input_aligned.pdb .
sbatch submit.sh
# Output: results/packed.pdb
```

Copy the output to `2dock/`:
```bash
cp results/packed.pdb ../2dock/packed.pdb
```

**For multi-target runs**, use `submit_all_prepack.py` from the GlycanDock root to submit all (GH, ligand) prepack jobs at once, automatically skipping any that already have `results/packed.pdb`:
```bash
python scripts/submit_all_prepack.py              # from GlycanDock root
python scripts/submit_all_prepack.py /path/to/GlycanDock
```

> ⚠️ **High `fa_rep` (~5000+ REU) in the prepack output is expected.** During prepack, the glycan is temporarily moved 1000 Å away and protein side chains repack into the empty binding site. When the glycan is placed back, interface clashes appear. GlycanDock Stage 2's soft-rep ramping resolves exactly these clashes — do not be alarmed by the energy.

---

## Step 4 — Docking Array (`4_dock/`)

**Run on cluster as a SLURM array. Each job takes ~30–60 minutes.**

```bash
# packed.pdb should already be here from Step 3
sbatch --array=1-200 submit.sh    # aim for ≥200 total structures
```

Each array task is an independent Monte Carlo trajectory producing one structure. `nstruct 1` per task is intentional — if a job is preempted on a `low` partition, only 1 structure is lost.

**Check the first log before submitting the full array:**
```
✓  Found 1 glycan trees          (not 2 or 3)
✓  No POLYMER_LOWER = 0 warning
✓  n_tor_cycles = 10
✓  n_tor_moves_accepted > 0
```

For multi-target runs, place `submit_dock.py` (from `scripts/`) in the GH folder and run it to submit all ligands at once — it automatically skips any that already have results:
```bash
cp scripts/submit_dock.py [GH]/submit_dock.py    # once per GH
python [GH]/submit_dock.py
```

---

## Step 5 — Score & Rank (`scripts/`)

**Run after all docking jobs complete — locally or on the cluster.**

```bash
# Option A — run locally
cd /path/to/your/GH/
python3 /path/to/scripts/glycandock_scores.py

# Option B — submit to cluster (Siegel lab: loads ligandmpnn_env conda environment)
cd /path/to/your/GH/
sbatch /path/to/scripts/Scoring.sh
```

The script scans all subdirectories for `results/` folders and applies a three-stage filter:

| Stage | Filter | Purpose |
|---|---|---|
| 1 | `n_tor_cycles ≥ 1` AND `ring_Lrmsd < 10 Å` | Remove invalid runs and structures that drifted out of the site |
| 2 | Top 20% by `total_score` | Keep low-energy structures |
| 3 | Top 10 by `interaction_energy` | Final ranking by binding energy |

Outputs are written to `Top_PDBs/` inside each docking subdirectory:
- `top_glycandock_summary.csv` — key scores for the top 10
- `top_glycandock_fullscores.csv` — all Rosetta score columns
- `top_glycandock_report.txt` — human-readable summary
- `*.pdb` — decompressed top 10 PDB structures

**Primary ranking metric: `interaction_energy`** (most negative = best binding)

---

## Score Columns Reference

| Column | Meaning |
|---|---|
| `interaction_energy` | Protein–glycan binding energy (REU); rank by this |
| `ring_Lrmsd` | Glycan displacement from input pose (Å); >10 Å = drifted out of site |
| `n_tor_cycles` | Torsion sampling cycles; 0 = wrong residue type (bad run) |
| `n_tor_moves_accepted` | Torsion moves accepted; should be > 0 |
| `total_score` | Overall Rosetta energy; used for pre-filtering only |
| `Fnat` | Fraction native contacts vs. input pose (only meaningful with crystal reference) |

---

## Input PDB Format Reference

The alignment scripts (Step 2) build this automatically. For reference:

```
HETNAM     Glc X   1  ->4)-beta-D-Glcp       ← tells Rosetta the correct residue type
HETNAM     Glc X   2  ->4)-beta-D-Glcp
HETNAM     Glc X   3  ->4)-beta-D-Glcp
LINK         O4  Glc X   1                 C1  Glc X   2                  1.50
LINK         O4  Glc X   2                 C1  Glc X   3                  1.50
ATOM      1  N   MET A   1  ...            ← protein, chain A
...
HETATM 1001  C1  Glc X   1  ...            ← glycan, chain X (MUST be at bottom)
...
TER
END
```

Key rules:
- Residue name: `Glc` (not BGC) — required with `-alternate_3_letter_codes pdb_sugar`
- LINK atoms: `O4` → `C1` (O4 is the glycosidic oxygen in Rosetta's Glc naming)
- Glycan HETATM records must come **after** all protein ATOM records
- One `TER` after the last HETATM, then `END`
- No CONECT records — they double-count with LINK records

---

## Troubleshooting

| Symptom | Root Cause | Fix |
|---|---|---|
| `sequence position requested was 0` + "Found 2 glycan trees" | `-auto_detect_glycan_connections` finds a spurious 3rd linkage by distance scan | Use `-maintain_links`; ensure LINK records are in PDB |
| `->3)-beta-D-Glcp` in log instead of `->4)` | Manual coordinate conversion distorted ring geometry | Use `pose_from_saccharide_sequence` (Step 1); never convert manually |
| `n_tor_cycles = 0` in score file | Wrong residue type assigned → torsion DOFs at C3 not C4 | Fix residue type first |
| Wrong residue type despite correct HETNAM | Missing HETNAM records or wrong format | Check HETNAM lines are present specifying `->4)-beta-D-Glcp` |
| Top structures show glycan outside active site | Missing `ring_Lrmsd` filter | Confirm `ring_Lrmsd < 10 Å` is in Stage 1 filter |
| Pair-fit RMSD > 0.5 Å | Atom mapping mismatch or docked pose has non-ideal ring geometry | Verify atom name correspondence; RMSD up to ~0.5 Å is acceptable |

---

## Saccharide Sequence Notation

Rosetta carbohydrate sequence string format: `linkage-monosaccharide-(bond->position)-...`

| Oligosaccharide | Abbreviation | Sequence string |
|---|---|---|
| Cellotriose (β-1→4 glucose ×3) | CL3 | `b-D-Glcp-(1->4)-b-D-Glcp-(1->4)-b-D-Glcp` |
| Laminaritrioside (β-1→3 glucose ×3) | CR3 | `b-D-Glcp-(1->3)-b-D-Glcp-(1->3)-b-D-Glcp` |
| Glc-Glc + α-Xyl-1,6 on reducing Glc | H2B | `[a-D-Xylp-(1->6)]b-D-Glcp-(1->4)-b-D-Glcp` |
| Glc-Glc-Glc + α-Xyl-1,6 on non-reducing Glc | H3B | `b-D-Glcp-(1->4)-b-D-Glcp-(1->4)-[a-D-Xylp-(1->6)]b-D-Glcp` |
| Xylotriose (β-1→4 xylose ×3) | XY3 | `b-D-Xylp-(1->4)-b-D-Xylp-(1->4)-b-D-Xylp` |

Abbreviations: `p` = pyranose; `a` = alpha, `b` = beta.

---

## Scripts Reference (`scripts/`)

| Script | Language | Purpose | Run from |
|---|---|---|---|
| `submit_all_prepack.py` | Python | Submit all `[GH]/[LIG]/1prepack/` jobs across the entire GlycanDock root; skips pairs with existing `results/packed.pdb` | GlycanDock root |
| `submit_dock.py` | Python | Submit all `[LIG]/2dock/` docking array jobs for one GH; skips ligands with existing results; checks `packed.pdb` exists | `[GH]/` folder (copy one here per GH) |
| `glycandock_scores.py` | Python | Load all `score_*.sc` files in `results/`, apply three-stage filter, write `Top_PDBs/` with CSVs, TXT report, and decompressed top 10 PDBs | `[GH]/` folder |
| `Scoring.sh` | SLURM bash | Cluster wrapper for `glycandock_scores.py`; loads conda `ligandmpnn_env` (Siegel lab path); submit with `sbatch Scoring.sh` | `[GH]/` folder |

**Three-stage filter in `glycandock_scores.py`:**

| Stage | Filter | Purpose |
|---|---|---|
| 1 | `n_tor_cycles ≥ 1` AND `ring_Lrmsd < 10 Å` | Remove invalid runs and structures that drifted out of the site |
| 2 | Top 20% by `total_score` | Keep low-energy structures |
| 3 | Top 10 by `interaction_energy` | Final ranking by binding energy |

---

## Citation

> Nance, M. L., et al. "GlycanDock: A Database-Guided Docking Algorithm for Oligosaccharide Ligands in Rosetta." *J. Chem. Theory Comput.* 17, 6799–6813 (2021). https://doi.org/10.1021/acs.jctc.1c00650
