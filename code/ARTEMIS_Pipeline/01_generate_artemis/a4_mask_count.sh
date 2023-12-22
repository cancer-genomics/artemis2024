#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_fsize=1000G
#$ -l mem_free=65G
#$ -l h_vmem=65G
#$ -l h_rt=96:00:00
#$ -o ./logs
#$ -pe local 4



jellyfish-linux count --mer-len 24 --size 10G --threads 4 --canonical --quality-start=32 --lower-count=1 --output masked.jellyfish masked.fasta

jellyfish-linux dump --column masked.jellyfish > masked.txt



