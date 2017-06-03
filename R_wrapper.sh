#!/bin/bash

dir=${1}
data_file=${1}${2}
analysis_file=${1}'analysis.txt'
output_file=${3}

source activate beests-hack
export LD_LIBRARY_PATH="/gpfs_home/dscott3/lib:"$LD_LIBRARY_PATH

Rscript run.R ${data_file} ${dir} ${analysis_file}
