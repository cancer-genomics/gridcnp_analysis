#!/bin/bash
#$ -cwd
#$ -l mem_free=10G
#$ -l h_vmem=10G

R CMD BATCH 3-normalizeToReference.R
