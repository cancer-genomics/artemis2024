#!/bin/bash
#$ -cwd
#$ -j y
#$ -R y
#$ -l mem_free=100G
#$ -l h_vmem=100G
#$ -l h_fsize=100G
#$ -l h_rt=24:00:00
#$ -o ./logs



rlib="${HOME}/Library/R/3.12-bioc-release-conda"
module load conda_R/4.0.x

CWD=$PWD
fragdir=$1

bindir=$2
binfile=$3
target=$4
mkdir -p $bindir

samplepath=$(find $fragdir -maxdepth 1 -name "*.rds"  | sort -u | head -n $SGE_TASK_ID | tail -n 1)
sample=$(basename $samplepath | awk '{ gsub(".rds", "") ; print $0}')

frag_file=${fragdir}/${sample}.rds
bin_file=${bindir}/${sample}_5mb.csv
if [ ! -f $bin_file ]; then
    R_LIBS_USER=$rlib Rscript x1_bin_corrected.r --id $frag_file --outdir $bindir --bins $binfile --target $target
fi ;
