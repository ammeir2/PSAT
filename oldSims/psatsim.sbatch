#!/bin/bash
#SBATCH --job-name fdrSim # Set a name for your job. This is useful if you have multiple jobs queued.
#SBATCH --partition short    # Slurm partition to use
#SBATCH --ntasks 1        # Number of tasks to run. By default, one CPU core will be allocated per task
#SBATCH --time 0-12:00       # Wall time limit in D-HH:MM
#SBATCH --mem-per-cpu=4000    # Memory limit for each tasks (in MB)
#SBATCH -o flowout.out     # File to which STDOUT will be written
#SBATCH -e flowout.out
#SBATCH --mail-type=END       # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=amitmeir@uw.edu # Email to which notifications will be sent
#SBATCH --array=1-10

module load R

R CMD BATCH --no-save --no-restore "--args setting='${SLURM_ARRAY_TASK_ID}'" mleSIM.R
