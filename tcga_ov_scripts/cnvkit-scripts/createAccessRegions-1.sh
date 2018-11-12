#!/bin/bash
#$ -cwd
#$ -j y
#$ -l mem_free=50G
#$ -l h_vmem=50G
#$ -l h_rt=48:00:00
#$ -l h_fsize=100G
#$ -l gwas

datadir=/dcl01/scharpf1/data/gridcnp_analysis/tcga_ov/cnvkit/data
mkdir -p $datadir
fasta=/dcl01/scharpf/data/reference/hg38/Homo_sapiens_assembly38.fasta
$cnvkit access $fasta -s 10000 -o $datadir/access-10kb.hg38.bed