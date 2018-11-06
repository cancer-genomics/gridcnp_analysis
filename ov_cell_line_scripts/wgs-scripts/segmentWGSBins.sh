#!/bin/bash
#$ -cwd
#$ -j y
#$ -l mem_free=10G
#$ -l h_vmem=10G
#$ -l h_rt=24:00:00
#$ -t 1-20

Rscript 1-segmentWGSBins.R $SGE_TASK_ID
