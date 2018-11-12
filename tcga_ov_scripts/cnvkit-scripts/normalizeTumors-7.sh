#!/bin/bash
#$ -cwd
#$ -j y
#$ -l mem_free=45G
#$ -l h_vmem=45G
#$ -l h_rt=48:00:00
#$ -l h_fsize=100G
#$ -t 1-11
#$ -l gwas

bamdir=/dcl01/scharpf1/data/bams/TCGA/download/OV
input=$(ls $bamdir | egrep '*-01[A-Z].bam$' | head -n $SGE_TASK_ID | tail -n 1)
bamfile=$bamdir/$input
sample=${input%.bam}

covdir=/dcl01/scharpf1/data/gridcnp_analysis/tcga_ov/cnvkit/tumorcov
reference=/dcl01/scharpf1/data/gridcnp_analysis/tcga_ov/cnvkit/reference/reference.cnn
outputdir=/dcl01/scharpf1/data/gridcnp_analysis/tcga_ov/cnvkit/tumornorm
mkdir -p $outputdir

$cnvkit fix \
    $covdir/$sample.targetcoverage.cnn $covdir/$sample.antitargetcoverage.cnn \
    $reference -o $outputdir/$sample.cnr
    
