#!/bin/bash
#SBATCH --job-name=shzBC
#SBATCH --time=04:00:00
#SBATCH --output=shzBC%A_%a.out
#SBATCH --error=shzBC%A_%a.err
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=32000
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

Rscript plot_curves.R
