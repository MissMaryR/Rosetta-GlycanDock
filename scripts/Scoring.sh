#!/bin/bash
#SBATCH --job-name=scoring
#SBATCH --output=logs/out_%A_%a.out 
#SBATCH --cpus-per-task=1
#SBATCH --partition=low
#SBATCH --time=1:00:00
#SBATCH --mem=5G

# Environment setup
export TORCH_HOME=/quobyte/jbsiegelgrp/software/LigandMPNN/.cache

module load conda/latest

eval "$(conda shell.bash hook)"

conda activate /quobyte/jbsiegelgrp/software/envs/ligandmpnn_env

# Run your Python scoring script
python /quobyte/jbsiegelgrp/missmaryr/scripts/glycandock_scores.py
