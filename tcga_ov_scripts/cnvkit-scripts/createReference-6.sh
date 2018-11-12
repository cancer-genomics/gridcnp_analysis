#!/bin/bash
#$ -cwd
#$ -j y
#$ -l mem_free=30G
#$ -l h_vmem=30G
#$ -l h_rt=48:00:00
#$ -l h_fsize=100G
#$ -l gwas

fasta=/dcl01/scharpf/data/reference/hg38/Homo_sapiens_assembly38.fasta
normals=/dcl01/scharpf1/data/gridcnp_analysis/tcga_ov/cnvkit/normalcov/*coverage.cnn
outdir=/dcl01/scharpf1/data/gridcnp_analysis/tcga_ov/cnvkit/reference
mkdir -p $outdir
$cnvkit reference $normals -f $fasta -o $outdir/reference.cnn