#!/bin/bash
#SBATCH --job-name=scoring
#SBATCH --output=logs/scoring_%j.out
#SBATCH --cpus-per-task=1
#SBATCH --partition=low
#SBATCH --time=1:00:00
#SBATCH --mem=5G

# Change to submission directory (where sbatch was called)
cd "$SLURM_SUBMIT_DIR"

# Create logs dir if it doesn't exist
mkdir -p logs

# Environment setup
export TORCH_HOME=/quobyte/jbsiegelgrp/software/LigandMPNN/.cache

module load conda/latest
eval "$(conda shell.bash hook)"
conda activate /quobyte/jbsiegelgrp/software/envs/ligandmpnn_env

# Score all ligands in this GH folder (finds */2dock/results/ automatically)
python /quobyte/jbsiegelgrp/missmaryr/scripts/glycandock_scores.py
