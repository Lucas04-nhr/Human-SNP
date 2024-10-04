#!/bin/bash
#SBATCH --job-name=preprocessing_GZ
#SBATCH --output=./log/Guangzhou/preprocessing_GZ.%j.out
#SBATCH --error=./log/Guangzhou/preprocessing_GZ.%j.err
#SBATCH --cpus-per-task=5
#SBATCH --mem=5G
#SBATCH --export=DATA_PATH='/mnt/raid6/bacphagenetwork/data/skin_metagenome/Guangzhou/02_rm_host',OUTPUT_BASE_PATH='/mnt/raid6/bacphagenetwork/data/'
#SBATCH --array=140-160%4

# Initialize the environment
echo "Initializing..."

# Set the paths of the output files

INDEXING_PATH="$OUTPUT_BASE_PATH/00_bwa_index/GRCh38"
INDEXING_FILE="$OUTPUT_BASE_PATH/00_bwa_index/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa"
ALIGNED_DATA_PATH="$OUTPUT_BASE_PATH/01_align/Guangzhou"
FILTERED_DATA_PATH="$OUTPUT_BASE_PATH/02_filter/Guangzhou"
SORTED_DATA_PATH="$OUTPUT_BASE_PATH/03_sort/Guangzhou"
CONVERTED_DATA_PATH="$OUTPUT_BASE_PATH/04_convert/Guangzhou"

echo "The original *.fastq files are located in $DATA_PATH."
echo "The indexing data is located in $INDEXING_PATH."
echo "The indexing genome data is $INDEXING_FILE."
echo "The alignment results will be saved in $ALIGNED_DATA_PATH."
echo "The *.bam files removed unmapped reads will be saved in $FILTERED_DATA_PATH."
echo "The sorted *.bam files and their indexes will be saved in $SORTED_DATA_PATH."
echo "The converted *.sam files will be saved in $CONVERTED_DATA_PATH."

echo "Initializing completed."
echo "=============================="

# Get the list of all *.fastq files
infile=($( cat GZ_sbatch.list | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print $0}' ))
file_base=$(basename "$infile" _1.fastq.gz)

# Extract the sample name
if [[ $infile == *Beijing* ]]; then
  sample_name=$(echo "$infile" | grep -oE 'BJ[0-9]{3}')
elif [[ $infile == *Guangzhou* ]]; then
  sample_name_temp=$(echo "$infile" | awk -F'/' '{n=split($NF,a,"_"); printf "%03d", a[n-1]}')
  sample_name="GZ$sample_name_temp"
else
  echo "Error: The sample name cannot be extracted."
  exit 1
fi

echo "Processing $sample_name..."

# Set the paths of the input files
GENE_DATA_1="$DATA_PATH/${file_base}_1.fastq.gz"
GENE_DATA_2="$DATA_PATH/${file_base}_2.fastq.gz"

# Step 1: Align the reads to the reference genome
echo "Aligning the reads to the reference genome..."
bwa mem -t 5 -R "@RG\tID:${sample_name}\tSM:${sample_name}\tPL:Illumina\tCN:GZ" $INDEXING_FILE $GENE_DATA_1 $GENE_DATA_2 -M > $ALIGNED_DATA_PATH/${sample_name}.sam \
|| { echo "Error: bwa mem failed in processing $sample_name."; exit 1; }
echo "Alignment completed."

# Step 2: Filter the unmapped reads
echo "Filtering the unmapped reads..."
samtools view -b -F 4 $ALIGNED_DATA_PATH/${sample_name}.sam > $FILTERED_DATA_PATH/${sample_name}.bam \
|| { echo "Error: samtools view failed in processing $sample_name."; exit 1; }
echo "Filtering completed."

# Step 3: Sort the filtered reads
echo "Sorting the filtered reads..."
samtools sort $FILTERED_DATA_PATH/${sample_name}.bam -o $SORTED_DATA_PATH/${sample_name}.bam \
|| { echo "Error: samtools sort failed in processing $sample_name."; exit 1; }
echo "Sorting completed."

echo "Indexing the sorted reads..."
samtools index $SORTED_DATA_PATH/${sample_name}.bam \
|| { echo "Error: samtools index failed in processing $sample_name."; exit 1; }
echo "Indexing completed."

# All done
echo "=============================="
echo "Processing $sample_name completed."
