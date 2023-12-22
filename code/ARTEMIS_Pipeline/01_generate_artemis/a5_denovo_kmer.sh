#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_fsize=1000G
#$ -l mem_free=65G
#$ -l h_vmem=65G
#$ -l h_rt=96:00:00
#$ -o ./logs
#$ -t 1-1287
#$ -pe local 4

name=$(find ./bed_type -maxdepth 1 -name "*.bed" | \
    sort -u | \
    head -n $SGE_TASK_ID | \
    tail -n 1)

name=$(basename $name)
sample=${name//.bed}
echo $sample

mkdir -p ./denovo_kmer_type

jellyfish-linux count --mer-len 24 --size 10G --threads 4 --canonical --quality-start=32 --lower-count=1 --output ./denovo_kmer_type/${sample}.jellyfish ./fasta_type/${sample}.fasta

jellyfish-linux dump --column ./denovo_kmer_type/${sample}.jellyfish > ./denovo_kmer_type/${sample}.txt


