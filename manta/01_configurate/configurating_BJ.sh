#!/bin/bash
#SBATCH --job-name=configurate_BJ
#SBATCH --output=./log/Beijing/configurate_BJ_%j.out
#SBATCH --error=./log/Beijing/configurate_BJ_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --export=INPUT_PATH='/mnt/raid6/bacphagenetwork/data/07_manta/00_format_converted/Beijing',OUTPUT_PATH='/mnt/raid6/bacphagenetwork/data/07_manta/01_exec/Beijing',MANTA_INSTALL_PATH='/mnt/raid6/bacphagenetwork/tools/manta/',REF_FILE_CH38='/mnt/raid6/bacphagenetwork/data/00_bwa_index/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa',REF_FILE_CHM13='/mnt/raid6/bacphagenetwork/data/00_bwa_index/chm13v2/chm13v2.0_noY.fa'
#SBATCH --array=2-201%5

# Initialize the environment
echo "Initializing..."

echo "The converted *.bam files and their indexes are located in $INPUT_PATH."
echo "The configuration files will be saved in $OUTPUT_PATH."

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
echo "Checking $INPUT_PATH/${sample_name}.bam..."
if [ ! -f "$INPUT_PATH/${sample_name}.bam" ]
then
    echo "Error: $INPUT_PATH/${sample_name}.bam does not exist."
    exit 2
else
    echo "The file $INPUT_PATH/${sample_name}.bam exists."
fi

echo "Checking $INPUT_PATH/${sample_nam}.bam.bai..."
if [ ! -f "$INPUT_PATH/${sample_name}.bam.bai" ]
then
    echo "Error: $INPUT_PATH/${sample_name}.bam.bai does not exist."
    exit 2
else
    echo "The file $INPUT_PATH/${sample_name}.bam.bai exists."
fi

# Check the reference files
echo "Checking the reference files..."
if [ ! -f "$REF_FILE_CH38" ]
then
    echo "Error: $REF_FILE_CH38 does not exist."
    exit 2
else
    echo "The reference file $REF_FILE_CH38 exists."
fi

if [ ! -f "$REF_FILE_CHM13" ]
then
    echo "Error: $REF_FILE_CHM13 does not exist."
    exit 2
else
    echo "The reference file $REF_FILE_CHM13 exists."
fi

echo "Initializing complete."
echo "=========================================================================="

# Configure the manta
echo "Configuring the manta..."
export MANTA_ANALYSIS_PATH=$OUTPUT_PATH/${sample_name}
mkdir -p $MANTA_ANALYSIS_PATH
echo "Now processing ${sample_name}..."
echo "The configuration files will be saved in $MANTA_ANALYSIS_PATH."

${MANTA_INSTALL_PATH}/bin/configManta.py \
--bam $INPUT_PATH/${sample_name}.bam \
--referenceFasta $REF_FILE_CHM13 \
--runDir ${MANTA_ANALYSIS_PATH} \
|| { echo "Error: configManta.py failed"; exit 3; }

echo "Configuring the manta complete."
echo "=========================================================================="

echo "The configuration files of ${sample_name} have been saved in $MANTA_ANALYSIS_PATH."
echo "All done."
