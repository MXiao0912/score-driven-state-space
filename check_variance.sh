#!/bin/bash
#SBATCH -D /rds/user/mx235/hpc-work/px_nav_ssm/mytvp/
#SBATCH -p sapphire
#SBATCH --job-name=verify_variance
#SBATCH --output=log/output_%a.out
#SBATCH --error=log/error_%A_%a.err
#SBATCH --array=1-100   # Define the range of tasks (0-9 means 10 tasks)
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

# Load necessary modules
module load matlab/r2024a

SAMPLE=$(head -n $SLURM_ARRAY_TASK_ID job.csv | tail -n 1 | tr -d '\n\r')

# Run the Python script with multiple arguments
srun matlab -nodisplay -nosplash -r "verify_variance('$SAMPLE'); exit force;"

# deactivate

