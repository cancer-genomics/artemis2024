#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_fsize=1000G
#$ -l mem_free=250G
#$ -l h_vmem=250G
#$ -l h_rt=96:00:00
#$ -o ./logs

module load conda_R
Rscript x1_split_bed.r
