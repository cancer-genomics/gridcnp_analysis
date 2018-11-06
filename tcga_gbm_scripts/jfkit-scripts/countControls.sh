#!/bin/bash

#$ -cwd
#$ -j y
#$ -R y
#$ -t 1-21
#$ -l mem_free=30G
#$ -l h_vmem=30G
#$ -l h_rt=96:00:00

WD=$PWD
bamdir=/dcl01/scharpf1/data/bams/TCGA/download/GBM
cd $bamdir
input=$(ls | egrep '*-10[A-Z].bam$' | head -n $SGE_TASK_ID | tail -n 1)
cd $WD
Rscript 2-countControls.R --id $input --bamdir $bamdir
