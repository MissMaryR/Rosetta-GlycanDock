#!/bin/bash --norc
#SBATCH --job-name=glycandock
#SBATCH --partition=low
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --output=logs/glycandock_%A_%a.out
#SBATCH --error=logs/glycandock_%A_%a.err
#SBATCH --requeue

# Step 4: GlycanDock Stage 1+2 docking array.
# Input:  packed.pdb  (prepacked output from 3_prepack/results/packed.pdb)
# Output: results/packed_<N>_0001.pdb.gz  (one per array task)
#         results/score_<N>.sc
#
# Submit: sbatch --array=1-200 submit.sh   (aim for >=200 total structures)
# nstruct 1 per array task is intentional — limits loss to 1 structure if preempted.
#
# After all jobs complete, run:  python3 ../script/glycandock_scores.py
#
# Edit ROSETTA_BIN and ROSETTA_DB for your cluster.

ROSETTA_BIN=/path/to/rosetta/main/source/bin
ROSETTA_DB=/path/to/rosetta/main/database

mkdir -p logs results

$ROSETTA_BIN/GlycanDock.static.linuxgccrelease \
  -database $ROSETTA_DB \
  @flags \
  -s packed.pdb \
  -nstruct 1 \
  -suffix _$SLURM_ARRAY_TASK_ID \
  -out:path:all results/ \
  -out:pdb_gz
