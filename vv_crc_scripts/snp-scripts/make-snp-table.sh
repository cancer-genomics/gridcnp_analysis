#! /bin/bash
#$ -cwd
#$ -l mem_free=20G
#$ -l h_vmem=20G
#$ -l h_rt=24:00:00
#$ -l h_fsize=100G
#$ -t 1-6

WD=$PWD

bamdir=/dcl01/scharpf1/data/bams/gridcnp_analysis/vv_crc/preprocessed_bams
cd $bamdir
tumor=$(ls -1v t*processed.bam | head -n $SGE_TASK_ID | tail -n 1)
normal=n${tumor:1}
cd $WD
output_dir=/dcl01/scharpf1/data/gridcnp_analysis/vv_crc/snps
mkdir -p $output_dir
sample=${tumor%.bam}
output_file=$output_dir/$sample.LOH.tsv

Rscript make-snp-table.R \
        -n $bamdir/$normal \
        -t $bamdir/$tumor \
        -c /dcl01/scharpf1/data/probe_sets/probe_bed_files/CpPa2_reduced_hg19.bed \
        -a hg19 \
        -o $output_file