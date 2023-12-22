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
#Where are the fastq
fastq=

#Where do you want the kmer counts to go
outdir=./counts

#Where are the bams
bams=

#where do you want the read counts to go
outdir_bam=./read_counts

#How many fastq are there
num_fastq=
#how many bams are there? It could be different 
num_bam=

###########################
#THESE ARE CHANGEABLE BUT PROBABLY DONT CHANGE THEM UNLESS YOU'RE DOING SOMETHING NONSTANDARD
#Where are all your kmers
ref_all=kmers_all.fasta

#Where are the type-specific kmers
ref_type=./final_kmers

##########################

#Getting the hashed kmers
qsub -l cancergen -t 1:$num_fastq a1_fast_count_samples.sh $fastq $outdir $ref_all

#Counting the read numbers for normalization -- if you don't have bams, don't do this
qsub -l cancergen -t 1:$num_bam a2_count_reads_denom.sh $bams $outdir_bam

#aggregate the kmer counts -- this will go to Hqw until fast_count is done
qsub -l cancergen -hold_jid a1_fast_count_samples.sh -t 1:$num_fastq a3_fast_aggregate.sh $outdir $ref_type

#Finalize
qsub -l cancergen -hold_jid a1_fast_count_samples.sh,a2_count_reads_denom.sh,a3_fast_aggregate.sh a4_feature_matrix.sh












