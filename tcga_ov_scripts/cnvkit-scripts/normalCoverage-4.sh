#!/bin/bash
#$ -cwd
#$ -j y
#$ -l mem_free=45G
#$ -l h_vmem=45G
#$ -l h_rt=48:00:00
#$ -l h_fsize=100G
#$ -t 1-20
#$ -l gwas

bamdir=/dcl01/scharpf1/data/bams/TCGA/download/OV
input=$(ls $bamdir | egrep '*-10[A-Z].bam$' | head -n $SGE_TASK_ID | tail -n 1)
bamfile=$bamdir/$input
sample=${input%.bam}

datadir=/dcl01/scharpf1/data/gridcnp_analysis/tcga_ov/cnvkit/normalcov
mkdir -p $datadir
targets=/dcl01/scharpf1/data/gridcnp_analysis/tcga_ov/cnvkit/targets/targets.bed
antitargets=/dcl01/scharpf1/data/gridcnp_analysis/tcga_ov/cnvkit/targets/antitargets.bed

$cnvkit coverage $bamfile $targets -o $datadir/$sample.targetcoverage.cnn
$cnvkit coverage $bamfile $antitargets -o $datadir/$sample.antitargetcoverage.cnn
