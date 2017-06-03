#!/bin/bash

### Necessary set-up ###
source activate beests-hack
export LD_LIBRARY_PATH="/gpfs_home/dscott3/lib:"$LD_LIBRARY_PATH

### Stuff for inspecting set-up ####
#conda list

#echo ' '
#echo $LD_LIBRARY_PATH | tr : '\n'
#echo ' '
#echo ${1}
#echo ${2}
#echo ${3}

### Actually function: MCMC SST Data ###
python run.py ${1} ${2} ${3}
Rscript run.R ${1} ${2} ${3}

