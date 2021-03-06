#!/bin/bash

#$ -cwd
#$ -j y
#$ -R y
#$ -t 1-6
#$ -l mem_free=10G
#$ -l h_vmem=10G
#$ -tc 65
#$ -l h_rt=96:00:00

WD=$PWD
bamdir=/dcl01/scharpf1/data/bams/gridcnp_analysis/vv_crc/preprocessed_bams
cd $bamdir
input=$(ls -1v t*processed.bam | head -n $SGE_TASK_ID | tail -n 1)
cd $WD
Rscript 1-countCases.R --id $input --bamdir $bamdir
