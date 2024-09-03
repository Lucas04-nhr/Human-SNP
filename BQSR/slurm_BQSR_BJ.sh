#!/bin/bash
#SBATCH --job-name=BQSR_BJ
#SBATCH --output=./log/01/Beijing/samtools_BJ_%j.out
#SBATCH --error=./log/01/Beijing/samtools_BJ_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=ANALYSIS_PATH='/mnt/raid6/bacphagenetwork/data/bwa_analysis/Beijing',INDEX_PATH='/mnt/raid6/bacphagenetwork/data/samtools_analysis/Beijing'
#SBATCH --array=1-201%4



gatk BaseRecalibrator -I dedup.bam -R reference.fasta --known-sites dbSNP.vcf.gz -O recal_data.table