#!/bin/bash --norc
#SBATCH --job-name=glycandock_prepack
#SBATCH --partition=low
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=logs/prepack_%j.out
#SBATCH --error=logs/prepack_%j.err

# Step 3: Stage 0 prepack — run ONCE before the docking array, takes ~2 min.
# Input:  input_aligned.pdb  (copy here from 0align/ output)
# Output: results/packed.pdb
#
# After this completes:
#   cp results/packed.pdb ../4_dock/packed.pdb
#   cd ../4_dock && sbatch submit.sh
#
# Edit ROSETTA_BIN and ROSETTA_DB for your cluster.

ROSETTA_BIN=/path/to/rosetta/main/source/bin
ROSETTA_DB=/path/to/rosetta/main/database

mkdir -p logs results

$ROSETTA_BIN/GlycanDock.static.linuxgccrelease \
  -database $ROSETTA_DB \
  @flags \
  -s input_aligned.pdb \
  -nstruct 1 \
  -out:path:all results/

# Rename output to packed.pdb
mv results/input_aligned_0001.pdb results/packed.pdb
