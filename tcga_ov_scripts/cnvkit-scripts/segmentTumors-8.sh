#!/bin/bash
#$ -cwd
#$ -j y
#$ -l mem_free=45G
#$ -l h_vmem=45G
#$ -l h_rt=48:00:00
#$ -l h_fsize=100G
#$ -t 1-11
#$ -l gwas

covdir=/dcl01/scharpf1/data/gridcnp_analysis/tcga_ov/cnvkit/tumornorm
input=$(ls $covdir | head -n $SGE_TASK_ID | tail -n 1)
covfile=$covdir/$input
sample=${input%.cnr}

outputdir=/dcl01/scharpf1/data/gridcnp_analysis/tcga_ov/cnvkit/tumorsegs
mkdir -p $outputdir

$cnvkit segment $covfile -o $outputdir/$sample.cns
    
