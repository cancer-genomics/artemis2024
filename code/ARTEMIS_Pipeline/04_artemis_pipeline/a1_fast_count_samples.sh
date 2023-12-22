#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_fsize=500G
#$ -l mem_free=26G
#$ -l h_vmem=26G
#$ -l h_rt=96:00:00
#$ -o ./logs
#$ -pe local 4

#---------------------------------------------------------------------
fastq=$1
outdir=$2
mkdir -p $outdir
refpath=$3

#--------------------------------------------------------------------------------

ref_name=$(basename $refpath | awk '{ gsub(".fasta", "") ; print $0 }')


samplepath=$(find $fastq -maxdepth 1 -name "*.fastq.gz" | \
    sed 's/\(.*\)_.*/\1/' | \
    sort -u | \
    head -n $SGE_TASK_ID | \
    tail -n 1)
sample=$(basename $samplepath | awk '{ gsub("_WGS", "") ; print $0 }')

fq1=$(find $fastq -maxdepth 1 -name "$sample*"|sort -u | head -n 1)
fq2=$(find $fastq -maxdepth 1 -name "$sample*"|sort -u | tail -n 1)

echo $sample
echo $outdir
echo $fastq

if [ -f ${outdir}/${sample}_${ref_name}_R1.jellyfish ]; then
   echo Sample $sample $ref_name has been processed. Exiting.
   exit 0
fi

echo $refpath
echo $ref_name
echo ${outdir}/${sample}_${ref_name}_R1.jellyfish
echo ${outdir}/${sample}/_${ref_name}_R2.jellyfish

jellyfish-linux count --mer-len 24 --size 16G --threads 4 --canonical --quality-start=32 --if $refpath --output ${outdir}/${sample}_${ref_name}_R1.jellyfish <(gunzip -c $fq1)

echo "one done!"
jellyfish-linux count --mer-len 24 --size 16G --threads 4 --canonical --quality-start=32 --if $refpath --output ${outdir}/${sample}_${ref_name}_R2.jellyfish  <(gunzip -c $fq2)

echo "two done!"
