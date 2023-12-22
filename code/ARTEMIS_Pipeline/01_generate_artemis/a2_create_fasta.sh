#!/bin/bash
#$ -cwd
#$ -j y
#$ -R y
#$ -l mem_free=50G
#$ -l h_vmem=50G
#$ -l h_fsize=100G
#$ -l h_rt=24:00:00
#$ -t 1-1287
#$ -o ./logs

module load bedtools
name=$(find ./bed_type -maxdepth 1 -name "*.bed" | \
    sort -u | \
    head -n $SGE_TASK_ID | \
    tail -n 1)
    
mkdir -p ./fasta_type

name=$(basename $name)
sample=${name//.bed}
bedtools getfasta  -fi chrm13.mod.fna -bed ./bed_type/$name -fo ./fasta_type/${sample}.fasta -name -fullHea