#!/bin/bash
#$ -cwd
#$ -j y
#$ -l mem_free=30G
#$ -l h_vmem=30G
#$ -l h_rt=48:00:00
#$ -l h_fsize=100G
#$ -l gwas

targets=/dcl01/scharpf1/data/bams/TCGA/capture/SureSelect_All_Exon_V2_with_annotation.hg38.bed
bam=/dcl01/scharpf1/data/bams/TCGA/download/OV/TCGA-04-1332-10A.bam
access=/dcl01/scharpf1/data/gridcnp_analysis/tcga_ov/cnvkit/data/access-10kb.hg38.bed
$cnvkit autobin $bam -t $targets -g $access