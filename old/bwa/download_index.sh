#! /bin/bash

echo "Downloading the genome data..."

# Set the download directory
echo "Please enter the path to the download directory:"
read download_dir
# Currently, the download directory is '/mnt/raid6/bacphagenetwork/data/bwa_index'
cd $download_dir

# wget -c https://processing.open-genomes.org/reference/CP086569.1-CHM13/CHM13_v1.1.fa

# wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_noY.fa.gz

wget -c http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz

echo "The genome data has been downloaded."

# Unzip the genome data
echo "Unzipping the genome data..."
gunzip chm13v2.0_noY.fa.gz
echo "The genome data has been unzipped."

# Print the path to the genome data
echo "The path to the genome data is $download_dir/chm13v2.0_noY.fa"