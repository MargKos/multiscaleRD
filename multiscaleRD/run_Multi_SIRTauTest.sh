#!/bin/bash
#SBATCH --job-name=MultiSIRTau                         # Job name
#SBATCH --partition=small                      # Partition name
#SBATCH --array=1-500                          
#SBATCH --mail-type=END                         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kostre@zib.de              # Where to send mail
#SBATCH --nodes=1                               # Run all processes on a single node
#SBATCH --ntasks=1                           # Number of tasks
#SBATCH --time=00-10:00:00                      # Time limit (necessary for Z1)
#SBATCH --output=job_%a-%a.log                  # Standard output and error log
#SBATCH --nodelist=htc-cmp[101-148]                    # Run on nodes between 101 and 148




date;hostname;pwd

cd /home/htc/bzfkostr/multiscaleRD/multiscaleRD
python3 Coupling_SIRTauTest.py $SLURM_ARRAY_TASK_ID



date