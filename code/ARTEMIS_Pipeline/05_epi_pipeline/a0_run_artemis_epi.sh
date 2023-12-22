#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_fsize=1000G
#$ -l mem_free=2G
#$ -l h_vmem=2G
#$ -l h_rt=96:00:00
#$ -o ./logs

module load conda_R

#############################
#YOU MUST SET THESE PATHS
#Where are the granges
granges= #Folder with the granges (this is an output of the DELFI pipeline in https://github.com/cancer-genomics/reproduce_lucas_wflow/tree/master/code/preprocessing)

#Where do you want the epi bins to go
outdir=./epi_bins

#How many samples
num_granges=386

#sequencing platform - can be hiseq "hi" or novaseq "nova", will depend on installation of
#https://github.com/cancer-genomics/PlasmaToolsNovaseq.hg19
#https://github.com/cancer-genomics/PlasmaToolsHiseq.hg19
#https://github.com/cancer-genomics/reproduce_lucas_wflow/tree/master/code/PlasmaTools.lucas
platform="hi"

###########################
#THESE ARE CHANGEABLE BUT PROBABLY DONT CHANGE THEM UNLESS YOU'RE DOING SOMETHING NONSTANDARD
#Where are all your kmers
ref_all=/dcs05/scharpf/data/annapragada/artemis_v2/finalize_references/Epi_Reference_Bins.csv

##########################

#Getting the bins
qsub -l cancergen -t 1:$num_granges a1_bin_corrected.sh $granges $outdir $ref_all $platform

#Finalize
qsub -l cancergen -hold_jid a1_bin_corrected.sh a2_matrix.sh












