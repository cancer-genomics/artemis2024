#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_fsize=1000G
#$ -l mem_free=10G
#$ -l h_vmem=10G
#$ -l h_rt=96:00:00
#$ -o ./logs

module load conda_R
Rscript x4_featureMatrix_unsc.r
