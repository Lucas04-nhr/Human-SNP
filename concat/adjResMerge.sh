#!/bin/bash

# 设置文件路径
input_directory="/mnt/raid6/bacphagenetwork/data/12_plink/Beijing/results/modified/c004"
output_file="/mnt/raid6/bacphagenetwork/data/12_plink/Beijing/results/modified/c004/merged_first_five_lines.csv"

# 初始化一个空的临时文件来存储合并的结果
temp_file=$(mktemp)

# 获取要处理的文件列表
input_files=$(ls $input_directory)

# 读取第一个文件的表头并写入结果文件
first_file=$(echo $input_files | awk '{print $1}')
head -n 1 "$input_directory/$first_file" > $temp_file

# 读取每个文件的前五行并合并
for input_file in $input_files; do
    tail -n +2 "$input_directory/$input_file" | head -n 5 >> $temp_file
done

# 将合并结果保存到新的文件
mv $temp_file $output_file