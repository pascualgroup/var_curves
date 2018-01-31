#!/bin/bash
#SBATCH --job-name=shzBC
#SBATCH --time=00:20:00
#SBATCH --output=slurm_output/shzBC_%A_%a.out
#SBATCH --error=slurm_output/shzBC_%A_%a.err
#SBATCH --array=1-50
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=3000
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

# 96%
Rscript accumulation_curves_midway.R BC 0.95 100 For_Shai_1099_isolate_42399_DBLa_types_at_96%_seq_ID.csv
# 90%
Rscript accumulation_curves_midway.R BC 0.95 100 For_Shai_1099_isolates_31840_types_at_90%_seq_ID.csv
# 80%
Rscript accumulation_curves_midway.R BC 0.95 100 For_Shai_1099_isolates_11546_types_at_80%_seq_ID.csv
# 70%
Rscript accumulation_curves_midway.R BC 0.95 100 For_Shai_1099_isolates_764_types_at_70%_seq_ID.csv
