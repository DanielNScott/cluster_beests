#!/bin/bash

# --------------------------------------------------------- #
# This file runs the BEESTS model and it's R post processing.
# --------------------------------------------------------- #

### Necessary set-up ###
# All-a-all-all-a-all beige erythang. By which I mean conda env.
source activate beests-hack

# Necessary for R to use appropriate libreadline.so
export LD_LIBRARY_PATH="/gpfs_home/dscott3/lib:"$LD_LIBRARY_PATH

### Stuff for inspecting set-up ####
#conda list

#echo ' '
#echo $LD_LIBRARY_PATH | tr : '\n'
#echo ' '
#echo ${1}
echo 'Running BEESTS model...'
echo 'Please note post_process.log will be overwritten.'

### Actually function: MCMC SST Data ###
python beests_model.py ${1}

if [ $? -ne 0 ]
then
   echo 'BEESTS terminated w/ non-zero exit status.'
   echo 'R post processing will not be run.'
	exit
fi

echo 'BEESTS exited with "success" status.'
echo 'Running R post processing.'
Rscript post_process.R ${1}
