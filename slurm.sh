#!/bin/bash

#SBATCH --job-name=nextstrain-a71   #This is the name of your job
#SBATCH --cpus-per-task=15                  #This is the number of cores reserved
#SBATCH --mem=30G              #Total memory
#Total memory reserved per core: 2GB

#SBATCH --time=06:00:00        #This is the time that your task will run: 06:00:00 
#SBATCH --qos=6hours           #You will run in this queue: 6hours

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=logs/slurm_output_%j.out      # Standard output (%j = job ID)
#SBATCH --error=logs/slurm_error_%j.err        # Standard error

#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=<your-mail>


#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.

# Load modules or activate conda/mamba environments
# eval "$(micromamba shell hook -s bash)"   # Ensure micromamba works inside batch jobs
# micromamba activate nextstrain            # Activate the Nextstrain environment

# Set environment variables
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export LMOD_DISABLE_SAME_NAME_AUTOSWAP="no"


# Run Snakemake
cd ~/nextstrain/enterovirus/enterovirus_a71/
time snakemake all_genes --cores $SLURM_CPUS_PER_TASK --nolock --stats run.stats --rerun-incomplete
