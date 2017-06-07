#!/bin/bash

###########################
## Slurm Config
###########################
#SBATCH -J beests-wtf
#SBATCH -t 02:00:00
#SBATCH -c 4
###########################

# Run the batch job.
dir=${1}
srun run_wrapper.sh ${dir} 1>${dir}/beests.out 2>>${dir}/beests.err

