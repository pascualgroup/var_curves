#!/bin/bash
#SBATCH --job-name=shzABC
#SBATCH --time=02:00:00
#SBATCH --output=shzABC%A_%a.out
#SBATCH --error=shzABC%A_%a.err
#SBATCH --array=1-20
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=8000
#SBATCH --partition=sandyb
##SBATCH --mail-type=END
##SBATCH --mail-user=pilosofs@uchicago.edu

# Print this sub-job's task ID
echo "Run ID: " $SLURM_ARRAY_TASK_ID

# Load modules
module load R
module load gcc/4.9

# general parameters
cd '/home/pilosofs/Shazia'

Rscript accumulation_curves_midway.R ABC $SLURM_ARRAY_TASK_ID 0.95
