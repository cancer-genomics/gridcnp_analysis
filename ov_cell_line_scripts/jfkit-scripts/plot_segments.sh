#!/usr/bin/bash

#$ -cwd
#$ -j y
#$ -l mem_free=2G
#$ -l h_vmem=2G
#$ -l h_rt=1:00:00
#$ -t 1-20

Rscript 4-plot_segments.R $SGE_TASK_ID