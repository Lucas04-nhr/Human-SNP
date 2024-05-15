#! /bin/bash

PWD=/mnt/raid6/bacphagenetwork/data/samtools_index

# Change the working directory
cd $PWD
echo "The working directory has been changed to $PWD."

# Download the reference genome
echo "Downloading the reference genome..."
wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz