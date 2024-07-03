#! /bin/bash
# Description: This script performs flagstat analysis on all BAM files in a specified directory and saves the results in a tab-separated text file.
# Set environment variables
export BAM_PATH='/mnt/raid6/bacphagenetwork/data/samtools_results/Guangzhou'
export OUTPUT_FILE='/mnt/raid6/bacphagenetwork/data/flagstat_summary/Guangzhou.txt'

# 检查samtools是否安装
if ! command -v samtools &> /dev/null
then
    echo "samtools could not be found, please install it first."
    exit 1
fi

# Check if the output directory exists
if [ ! -d $(dirname $OUTPUT_FILE) ]; then
    echo "Output directory $(dirname $OUTPUT_FILE) does not exist. Creating..."
    mkdir -p $(dirname $OUTPUT_FILE)
fi

# Check if the output file already exists
if [ -f $OUTPUT_FILE ]; then
    echo "Error: Output file $OUTPUT_FILE already exists. Would you like to overwrite it? (y/n)"
    read response
    if [ "$response" != "y" ]; then
        echo "Exiting script."
        exit 1
    fi
fi

# 初始化输出文件标题行
echo -e "Sample\tTotal\tSecondary\tSupplementary\tDuplicates\tMapped\tPaired in sequencing\tRead1\tRead2\tProperly paired\tWith itself and mate mapped\tSingletons\tWith mate mapped to a different chr\tWith mate mapped to a different chr (mapQ>=5)" > $OUTPUT_FILE

# 遍历BAM_PATH下的所有BAM文件
for bamfile in $BAM_PATH/*.bam
do
    # 获取文件名作为样本名
    sample_name=$(basename $bamfile .bam)

    # Print the sample name
    echo "Processing $sample_name..."

    # 运行samtools flagstat并提取信息
    stats=$(samtools flagstat $bamfile)
    
    # 解析flagstat输出并格式化
    total=$(echo "$stats" | grep "in total" | cut -d ' ' -f 1)
    secondary=$(echo "$stats" | grep "secondary" | cut -d ' ' -f 1)
    supplementary=$(echo "$stats" | grep "supplementary" | cut -d ' ' -f 1)
    duplicates=$(echo "$stats" | grep "duplicates" | cut -d ' ' -f 1)
    mapped=$(echo "$stats" | grep "mapped (" | cut -d ' ' -f 1)
    paired=$(echo "$stats" | grep "paired in sequencing" | cut -d ' ' -f 1)
    read1=$(echo "$stats" | grep "read1" | cut -d ' ' -f 1)
    read2=$(echo "$stats" | grep "read2" | cut -d ' ' -f 1)
    properly_paired=$(echo "$stats" | grep "properly paired" | cut -d ' ' -f 1)
    with_itself=$(echo "$stats" | grep "with itself and mate mapped" | cut -d ' ' -f 1)
    singletons=$(echo "$stats" | grep "singletons (" | cut -d ' ' -f 1)
    diff_chr=$(echo "$stats" | grep "with mate mapped to a different chr" | cut -d ' ' -f 1)
    diff_chr_mapq5=$(echo "$stats" | grep "with mate mapped to a different chr (mapQ>=5)" | cut -d ' ' -f 1)
    
    # 将结果写入输出文件
    echo -e "$sample_name\t$total\t$secondary\t$supplementary\t$duplicates\t$mapped\t$paired\t$read1\t$read2\t$properly_paired\t$with_itself\t$singletons\t$diff_chr\t$diff_chr_mapq5" >> $OUTPUT_FILE
done

echo "Flagstat analysis completed. Results are saved in $OUTPUT_FILE."