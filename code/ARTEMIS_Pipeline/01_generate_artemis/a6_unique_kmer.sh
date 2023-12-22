#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_fsize=200G
#$ -l mem_free=250G
#$ -l h_vmem=250G
#$ -l h_rt=96:00:00
#$ -o ./logs
#$ -t 1-1287
module load conda_R
mkdir -p ./unique_kmer_type

Rscript x6_select_kmers.r $SGE_TASK_ID

