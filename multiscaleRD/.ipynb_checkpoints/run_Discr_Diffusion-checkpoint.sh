#!/bin/bash
#SBATCH --job-name=DiscDiffE                        # Job name
#SBATCH --partition=small                      # Partition name                          
#SBATCH --mail-type=END                         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kostre@zib.de              # Where to send mail
#SBATCH --nodes=1                               # Run all processes on a single node
#SBATCH --ntasks=1                           # Number of tasks
#SBATCH --time=00-02:00:00                      # Time limit (necessary for Z1)
#SBATCH --output=job_%a-%a.log                  # Standard output and error log

date;hostname;pwd

cd /home/htc/bzfkostr/multiscaleRD/multiscaleRD
python3 Discretization_DiffusionExplicit.py


date
