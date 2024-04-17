#! /bin/bash

echo "Downloading the genome data..."

# Set the download directory
echo "Please enter the path to the download directory:"
read download_dir
# Currently, the download directory is '/mnt/raid6/bacphagenetwork/data/bwa_index'
cd $download_dir

wget -c https://public-ehpc-package.oss-cn-hangzhou.aliyuncs.com/lifescience/b37_human_g1k_v37.fasta
wget -c https://public-ehpc-package.oss-cn-hangzhou.aliyuncs.com/lifescience/b37_1000G_phase1.indels.b37.vcf
wget -c https://public-ehpc-package.oss-cn-hangzhou.aliyuncs.com/lifescience/b37_Mills_and_1000G_gold_standard.indels.b37.vcf
wget -c https://public-ehpc-package.oss-cn-hangzhou.aliyuncs.com/lifescience/b37_dbsnp_138.b37.vcf.zip
unzip b37_dbsnp_138.b37.vcf.zip

wget -c https://public-ehpc-package.oss-cn-hangzhou.aliyuncs.com/lifescience/gatk-examples_example1_NA20318_ERR250971_1.filt.fastq.gz
gunzip gatk-examples_example1_NA20318_ERR250971_1.filt.fastq.gz

wget -c https://public-ehpc-package.oss-cn-hangzhou.aliyuncs.com/lifescience/gatk-examples_example1_NA20318_ERR250971_2.filt.fastq.gz
gunzip gatk-examples_example1_NA20318_ERR250971_2.filt.fastq.gz

echo "The genome data has been downloaded."