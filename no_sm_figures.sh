#!/bin/sh

# make sure it's to Lar's account
#SBATCH --account larsms
# set time
#SBATCH --time=1-03:00:00

#set a name for the job visible in 'squeue'
#SBATCH --job-name="snakejob"

# set nodes
#SBATCH --nodes=1

# set one task
#SBATCH --ntasks=7

# one cpu/core per task
#SBATCH --cpus-per-task=1

#OPTIONAL
# who to send mail to
#SBATCH --mail-user=kwells4@stanford.edu
# what type of mail to send
#SBATCH --mail-type=END

. /home/kwells4/anaconda3/etc/profile.d/conda.sh
conda activate /labs/larsms/projects/mTEC_dev/mtec_snakemake/.snakemake/conda/1bdf5e39

Rscript scripts/figures_no_sm.R