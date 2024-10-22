#!/bin/bash
#SBATCH -D /rds/user/mx235/hpc-work/px_nav_ssm/mytvp/mymodel/tvp/
#SBATCH -p icelake
#SBATCH --job-name=verify_variance
#SBATCH --output=log/output_%a.out
#SBATCH --error=log/error_%A_%a.err
#SBATCH --array=1-137   # Define the range of tasks (0-9 means 10 tasks) 1-137
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# Load necessary modules
module load matlab/r2024a

# Run the script with multiple arguments
srun matlab -nodisplay -nosplash -r "rerun_diagnostics('$SLURM_ARRAY_TASK_ID'); exit force;"

# deactivate
