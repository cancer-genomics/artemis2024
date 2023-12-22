#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_fsize=1000G
#$ -l mem_free=200G
#$ -l h_vmem=200G
#$ -l h_rt=96:00:00
#$ -t 1-1287
#$ -o ./logs

module load conda_R
mkdir -p ./kmers_non_masked
Rscript x7_select_nonmask_kmers.r $SGE_TASK_ID
