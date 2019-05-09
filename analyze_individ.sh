ure it's to Lar's account
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

snakemake -j 999 --use-conda --cluster-config cluster.json --cluster "sbatch --account {cluster.account} --partition {cluster.partition} --job-name {cluster.job-name} --time {cluster.time} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus-per-task} --nodes {cluster.nodes} --mem {cluster.mem} --mail-user {cluster.mail} --mail-type {cluster.type}"

