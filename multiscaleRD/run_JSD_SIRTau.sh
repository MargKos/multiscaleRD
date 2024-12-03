<#!/bin/bash
#SBATCH --job-name=JSDSIRTau                        # Job name
#SBATCH --partition=small                       # Partition name
#SBATCH --mail-type=END                         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kostre@zib.de               # Where to send mail
#SBATCH --nodes=1                               # Run all processes on a single node
#SBATCH --ntasks=1                              # Number of tasks
#SBATCH --time=00-10:00:00                      # Time limit (days-hours:minutes:seconds)
#SBATCH --output=job_%j.log                     # Standard output and error log (%j is the job ID)

# Print the current date, hostname, and working directory
date
hostname
pwd

# Change to the specified directory
cd /home/htc/bzfkostr/multiscaleRD/multiscaleRD

# Run the Python script

python3 JSD_Run_SIRTauFast.py $SLURM_ARRAY_TASK_ID




