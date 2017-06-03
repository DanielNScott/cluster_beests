#!/bin/bash

#SBATCH -J beests-wtf
#SBATCH -t 02:00:00
#SBATCH -c 4

dir=${1}
data_file=${1}${2}
analysis_file=${1}'analysis.txt'
output_file=${3}

srun run_wrapper.sh ${data_file} ${dir} ${analysis_file} 1>>${output_file} 2>>${output_file}

