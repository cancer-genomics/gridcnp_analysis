#!/bin/bash

#$ -cwd
#$ -l mem_free=5G
#$ -l h_vmem=5G

R CMD BATCH 2-normalizeToReference.R
