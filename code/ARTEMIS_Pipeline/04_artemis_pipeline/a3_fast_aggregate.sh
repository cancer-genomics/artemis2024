#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_fsize=500G
#$ -l mem_free=20G
#$ -l h_vmem=20G
#$ -l h_rt=96:00:00
#$ -o ./logs

#---------------------------------------------------------------------
#The two paths below are the same, replace for the path to your kmer counts. -t above is the number of samples
path=$1
outdir=$1
ref=$2

#--------------------------------------------------------------------------------
samplepath=$(find $path -maxdepth 1 -name "*.jellyfish" | \
    sed 's/\_kmers_all.*//' | \
    sort -u | \
    head -n $SGE_TASK_ID | \
    tail -n 1)
sample=$(basename $samplepath | awk '{ gsub("_WGS", "") ; print $0 }')

echo $sample
echo $outdir

mkdir -p $outdir/features

if [ -f ${outdir}/features/${sample}.txt ]; then
   echo Sample $sample_R1 has been processed. Exiting.
   exit 0
fi

echo $path/${sample}_kmers_all_R1.jellyfish
echo $path/${sample}_kmers_all_R2.jellyfish

for i in {1..1280}
do
	refpath=$(find $ref -maxdepth 1 -name "*.fasta" | \
    sort -u | \
    head -n $i | \
    tail -n 1)
	ref_name=$(basename $refpath | awk '{ gsub(".fasta", "") ; print $0 }')
	echo $ref_name
	r1=$(jellyfish-linux query $path/${sample}_kmers_all_R1.jellyfish -s $refpath | awk '{ sum += $2 } END { print sum }')
	r2=$(jellyfish-linux query $path/${sample}_kmers_all_R2.jellyfish -s $refpath | awk '{ sum += $2 } END { print sum }')
	val=$(($r1 + $r2))
	echo $ref_name $val $sample >> ${outdir}/features/$sample.txt
	echo $i
done
