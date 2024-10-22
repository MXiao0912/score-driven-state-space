#!/bin/bash
#SBATCH -D /rds/user/mx235/hpc-work/px_nav_ssm/mytvp/mymodel/
#SBATCH -p icelake
#SBATCH --job-name=verify_variance
#SBATCH --output=log/output_%a.out
#SBATCH --error=log/error_%A_%a.err
#SBATCH --array=1-137   # Define the range of tasks (0-9 means 10 tasks)
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# Load necessary modules
module load matlab/r2024a

# Run the script with multiple arguments
srun matlab -nodisplay -nosplash -r "run_tvp_mymodel_tot('$SLURM_ARRAY_TASK_ID'); exit force;"

# deactivate
