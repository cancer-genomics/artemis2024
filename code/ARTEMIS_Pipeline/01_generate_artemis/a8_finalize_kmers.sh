#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_fsize=1000G
#$ -l mem_free=250G
#$ -l h_vmem=250G
#$ -l h_rt=96:00:00
#$ -t 1-1287
#$ -o ./logs

module load conda_R
mkdir -p ./final_kmers
Rscript x8_finalize_kmer_set.r $SGE_TASK_ID
