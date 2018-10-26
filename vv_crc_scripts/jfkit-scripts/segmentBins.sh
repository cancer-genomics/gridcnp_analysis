#!/bin/bash

#$ -cwd
#$ -j y
#$ -l mem_free=3G
#$ -l h_vmem=3G
#$ -l h_rt=1:00:00
#$ -t 1-6

Rscript 4-segmentBins.R $SGE_TASK_ID
