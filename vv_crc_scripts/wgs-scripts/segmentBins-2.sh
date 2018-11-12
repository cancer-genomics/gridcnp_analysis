#!/bin/bash
#$ -cwd
#$ -j y
#$ -R y
#$ -t 1-6
#$ -l mem_free=10G
#$ -l h_vmem=10G
#$ -l h_rt=48:00:00
Rscript 2-segmentBins.R $SGE_TASK_ID