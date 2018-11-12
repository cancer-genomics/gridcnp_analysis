#!/bin/bash
#$ -cwd
#$ -j y
#$ -R y
#$ -t 1-11
#$ -l mem_free=40G
#$ -l h_vmem=40G
#$ -l h_rt=96:00:00

WD=$PWD
bamdir=/dcl01/scharpf1/data/bams/TCGA/download/OV
cd $bamdir
input=$(ls | egrep '*-01[A-Z].bam$' | head -n $SGE_TASK_ID | tail -n 1)
cd $WD
Rscript 1-countCases.R --id $input --bamdir $bamdir
