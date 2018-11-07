#!/bin/bash
#$ -cwd
#$ -j y
#$ -R y
#$ -t 1-6
#$ -l mem_free=70G
#$ -l h_vmem=70G
#$ -l h_rt=48:00:00

WD=$PWD
# Change lines below after mark duplicates has finished running
#bamdir=/dcl01/scharpf1/data/bams/colorectal/tissue/wholegenome/preprocessed_bams
bamdir=/dcl01/scharpf1/data/bams/colorectal/tissue/wholegenome
cd $bamdir
#sample=$(ls -1v *processed.bam | head -n $SGE_TASK_ID | tail -n 1)
sample=$(ls -1v *sortednofa.bam | head -n $SGE_TASK_ID | tail -n 1)
input=$bamdir/$sample
cd $WD
Rscript 1-wgsReadDepth.R $input
