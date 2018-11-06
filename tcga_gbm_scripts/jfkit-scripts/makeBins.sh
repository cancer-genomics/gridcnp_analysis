#!/bin/bash

#$ -cwd
#$ -j y
#$ -l mem_free=10G
#$ -l h_vmem=10G

R CMD BATCH 0-makeBins.R
