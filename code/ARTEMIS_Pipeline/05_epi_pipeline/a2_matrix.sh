#!/bin/bash
#$ -cwd
#$ -j y
#$ -R y
#$ -l mem_free=10G
#$ -l h_vmem=10G
#$ -l h_fsize=10G
#$ -l h_rt=24:00:00
#$ -o ./logs


module load conda_R/4.0.x

Rscript x2_aggregate_matrix.r
