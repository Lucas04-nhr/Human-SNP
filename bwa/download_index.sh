#! /bin/bash

echo "Downloading the genome data..."

# Set the download directory
echo "Please enter the path to the download directory:"
read download_dir
# Currently, the download directory is '/mnt/raid6/bacphagenetwork/data/bwa_index'
cd $download_dir

wget -c https://processing.open-genomes.org/reference/CP086569.1-CHM13/CHM13_v1.1.fa

echo "The genome data has been downloaded."