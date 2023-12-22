#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_fsize=5G
#$ -l mem_free=2G
#$ -l h_vmem=2G
#$ -l h_rt=144:00:00
#$ -t 1-226
#$ -o ./logs

module load conda_R/4.1.x

#fold="fold"$SGE_TASK_ID

Rscript Ensemble_ARTEMIS1_delfi_LOO.r $SGE_TASK_ID