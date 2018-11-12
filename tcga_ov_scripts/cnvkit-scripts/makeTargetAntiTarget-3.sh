#!/bin/bash
#$ -cwd
#$ -j y
#$ -l mem_free=30G
#$ -l h_vmem=30G
#$ -l h_rt=48:00:00
#$ -l h_fsize=100G
#$ -l gwas

datadir=/dcl01/scharpf1/data/gridcnp_analysis/tcga_ov/cnvkit/targets
mkdir -p $datadir
targets=/dcl01/scharpf1/data/bams/TCGA/capture/SureSelect_All_Exon_V2_with_annotation.hg38.bed
access=/dcl01/scharpf1/data/gridcnp_analysis/tcga_ov/cnvkit/data/access-10kb.hg38.bed
$cnvkit target $targets --split -o $datadir/targets.bed
$cnvkit antitarget $targets -g $access -a 50000 -o $datadir/antitargets.bed
