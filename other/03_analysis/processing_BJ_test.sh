#!/bin/bash
#SBATCH --job-name=analysis_BJ
#SBATCH --output=./log/Beijing/analysis_BJ_%j.out
#SBATCH --error=./log/Beijing/analysis_BJ_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=INPUT_PATH='/mnt/raid6/bacphagenetwork/data/06_unmapped_removed/Beijing',OUTPUT_PATH='/mnt/raid6/bacphagenetwork/data/08_analysis/Beijing'
#SBATCH --array=1

# Initialate conda
source ~/.bashrc

conda activate snp_analysis

echo "The *.sam files whose unmapped area has been removed are located in $INPUT_PATH."
echo "The analysis files will be saved in $OUTPUT_PATH."

infile=($( cat BJ_sbatch.list | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print $0}' ))
sample_name=$(echo "$infile" | grep -oE 'BJ[0-9]{3}')

# Check all the needed files exist
# Check the output folder
if [ ! -d "$OUTPUT_PATH" ]
then
    echo "The output folder does not exist."
    exit 1
fi

# Check the input file
echo "Checking $INPUT_PATH/${sample_name}.removed.sam..."
if [ ! -f "$INPUT_PATH/${sample_name}.removed.sam" ]
then
    echo "Error: $INPUT_PATH/${sample_name}.removed.sam does not exist."
    exit 2
else
    echo "The file $INPUT_PATH/${sample_name}.removed.sam exists."
fi

# Analysis
echo "Analysis for $INPUT_PATH/${sample_name}.removed.sam..."
echo "========================================================================"

log_file=/mnt/raid6/bacphagenetwork/niehaoran/Human-SNP/other/03_analysis/log/Beijing/analysis_BJ_${SLURM_JOB_ID}.out
python /mnt/raid6/bacphagenetwork/niehaoran/Human-SNP/other/03_analysis/processing.py \
    -i $INPUT_PATH/${sample_name}.removed.sam \
    -o $OUTPUT_PATH/ \
    -s $OUTPUT_PATH/static/ \
>> $log_file 2>&1 \
|| { echo "Error: processing of ${sample_name} failed"; exit 4; }

echo "========================================================================"

echo "The analysis for $INPUT_PATH/${sample_name}.removed.sam is done."

# If SLURM reaches the last sample, delete the tmp folder
if [ ${SLURM_ARRAY_TASK_ID} -eq 1 ]
then
    rm -rf $OUTPUT_PATH/tmp
fi
