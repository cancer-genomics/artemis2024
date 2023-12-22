#!/bin/bash
#$ -cwd
#$ -j y
#$ -R y
#$ -l mem_free=50G
#$ -l h_vmem=50G
#$ -l h_fsize=100G
#$ -l h_rt=24:00:00
#$ -o ./logs

module load bedtools

#This is to create a version of chm13 with all of the repeat masked sections masked

bedtools maskfasta -fi chrm13.mod.fna -bed ucsc-t2t-repeat-masker.bed -fo masked.fasta
