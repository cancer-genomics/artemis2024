#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_fsize=100G
#$ -l mem_free=25G
#$ -l h_vmem=25G
#$ -l h_rt=96:00:00
#$ -o ./logs


module load samtools

# Inputs
dir=$1
outdir=$2
mkdir -p $outdir

# Main
samplepath=$(find $dir -maxdepth 1 -name "*.bam" | \
    sort -u | \
    head -n $SGE_TASK_ID | \
    tail -n 1)
sample=$(basename $samplepath | awk '{ gsub(".sorted", "") ; print $0 }')

sample=${sample//.bam}

if [ ! -f $outdir/$sample.txt ]
then
	echo $sample > $outdir/$sample.txt
	echo $sample
	samtools view -c -q 30 -F 3844 $samplepath >> $outdir/$sample.txt
else
	echo DONE!
fi

