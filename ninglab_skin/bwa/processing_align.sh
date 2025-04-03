#! /bin/bash
#SBATCH --job-name=ninglab_data_align
#SBATCH --output=./ninglab_data_align.%j.out
#SBATCH --error=./ninglab_data_align.%j.err
#SBATCH --cpus-per-task=4
#SBATCH --array=1-201%4
#SBATCH --mem=64G

# Initialize the environment
echo "Initializing the environment..."
echo "=============================="

# Load conda
echo "Activating conda environment..."
source /home/bacphagenetwork/.bashrc
conda activate snp_analysis
echo "The conda environment has been activated."
echo "=============================="


# Set the path to the genome data
GENOME_PATH='/mnt/raid6/bacphagenetwork/data/ninglab_skin'
INDEXING_PATH='/mnt/raid6/bacphagenetwork/data/00_bwa_index/GRCh38'
INDEXING_FILE='/mnt/raid6/bacphagenetwork/data/00_bwa_index/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa'
ANALYSIS_PATH='/mnt/raid6/bacphagenetwork/data/ninglab/01_bwa_analysis'

# Prompt the paths
echo "The path to the genome data has been set to $GENOME_PATH."
echo "The indexing file is located at $INDEXING_FILE."
echo "The analysis path is set to $ANALYSIS_PATH."
echo "Checking whether the analysis path exists..."
  
if [ -d "$ANALYSIS_PATH" ]; then
    echo "The analysis path exists."
else
    echo "The analysis path does not exist. Creating"
    mkdir -p "$ANALYSIS_PATH"
    echo "The analysis path has been created."
fi
echo "=============================="

# Analysing the data
echo "Analysing the genome data..."
infile=($( cat sbatch.list | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print $0}' ))
sample_name=$(basename "${infile}" | awk -F'_' '{print $1}')
echo "Processing ${sample_name}..."
file_fq1="$GENOME_PATH/${sample_name}_R1.fq.gz"
file_fq2="${GENOME_PATH}/${sample_name}_R2.fq.gz"
echo -e "Files for ${sample_name} are:\n ${file_fq1} and\n ${file_fq2}."
touch $ANALYSIS_PATH/${sample_name}.sam
mkdir -p $ANALYSIS_PATH/log
bwa mem -t 4 $INDEXING_PATH $file_fq1 $file_fq2 -a > $ANALYSIS_PATH/log/${sample_name}.log 2>&1  > $ANALYSIS_PATH/${sample_name}.sam || { echo "Error: bwa mem failed in processing $sample_name."; exit 1; }
echo "The analysis of $sample_name has been completed."
