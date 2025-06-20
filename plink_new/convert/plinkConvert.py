import pandas as pd
import numpy as np
import argparse
import os

# 设置命令行参数解析
parser = argparse.ArgumentParser(description='Merge and filter files in a directory.')
parser.add_argument('-i', '--input-directory', type=str, required=True, help='The directory containing the files to be processed')
parser.add_argument('-o', '--output-directory', type=str, required=True, help='The output directory for the filtered results')
parser.add_argument('-p', '--pheno-file', type=str, required=True, help='The phenotype file path for the analysis')
args = parser.parse_args()

# 获取目录路径和输出文件路径
input_directory = args.input_directory # /mnt/raid6/bacphagenetwork/data/12_plink_Full/results
output_directory = args.output_directory # /mnt/raid6/bacphagenetwork/data/12_plink_Full/converted
pheno_file = args.pheno_file # /mnt/raid6/bacphagenetwork/data/12_plink_Full/phenotype_full.tsv

if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# 读取表型文件
try:
    pheno_df = pd.read_csv(pheno_file, sep='\t')
    print(f"Phenotype file loaded with {len(pheno_df)} rows and {len(pheno_df.columns)} columns.")
except FileNotFoundError as fnfe:
    raise FileNotFoundError(f"Error while reading the phenotype file: {fnfe}")

count = 0
# 遍历输入目录下的所有 *.adjusted 文件
for filename in os.listdir(input_directory):
    count += 1
    if filename.endswith('.adjusted'):
        pheno_number = filename.split('.')[1][1:]  # 提取文件名中的数字部分
        file_path = os.path.join(input_directory, filename)
        df = pd.read_csv(file_path, sep='\s+')
        print(f"Processing file {filename} with {len(df)} rows and {len(df.columns)} columns.")
        # 增加"Bacterium"列
        df['Bacterium'] = pheno_df.columns[int(pheno_number) + 1] # index begins from 0
        df.to_csv(os.path.join(output_directory, filename), index=False)
        print(f"File {filename} processed and saved to {output_directory}.")
        print(f"Processed {count} files so far.")

print(f"All files processed. Total files processed: {count}.")
print(f"Filtered files saved to {output_directory}.")
