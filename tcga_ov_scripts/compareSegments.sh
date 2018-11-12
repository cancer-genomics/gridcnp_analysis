#!/usr/bin/bash

#$ -cwd
#$ -j y
#$ -l mem_free=5G
#$ -l h_vmem=5G
#$ -l h_rt=1:00:00
#$ -t 1-11

Rscript compareSegments.R $SGE_TASK_ID