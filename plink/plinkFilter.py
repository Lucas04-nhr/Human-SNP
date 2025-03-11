import os
import pandas as pd
import argparse
import sys

# 设置命令行参数解析
parser = argparse.ArgumentParser(description='Merge and filter files in a directory.')
parser.add_argument('-d', '--directory', type=str, required=True, help='The directory containing the files to be processed')
parser.add_argument('-o', '--output_file', type=str, required=True, help='The output file path for the filtered results')
parser.add_argument('-i', '--item', type=str, required=True, help='The item to filter on')
parser.add_argument('-t', '--threshold', type=float, required=True, help='The threshold value for filtering')
args = parser.parse_args()

# 获取目录路径和输出文件路径
directory = args.directory
output_file = args.output_file
threshold = args.threshold

# 设置可排序的列名
sortable_columns = ['UNADJ', 'GC', 'BONF', 'HOLM', 'SIDAK_SS', 'SIDAK_SD', 'FDR_BH', 'FDR_BY']

if args.item not in sortable_columns:
    raise ValueError(f"Invalid item to filter on: {args.item}. Please choose from {sortable_columns}")
else:
    sort_columns = args.item

# 初始化一个空的DataFrame用于存储合并后的数据
combined_df = pd.DataFrame()

try:
    # 遍历目录下的所有文件
    for filename in os.listdir(directory):
        if filename.endswith('.adjusted'):
            file_path = os.path.join(directory, filename)
            df = pd.read_csv(file_path, sep='\s+')
            print(f"Processing file {filename} with {len(df)} rows and {len(df.columns)} columns.")
            combined_df = pd.concat([combined_df, df], ignore_index=True, sort=False)
except FileNotFoundError as fnfe:
    raise FileNotFoundError(f"Error while processing files: {fnfe}")

try:
    # 筛选出FDR_BY列小于5 * 10e-8的行
    filtered_df = combined_df[combined_df[sort_columns] <= threshold]
    print(f"Filtered data contains {len(filtered_df)} rows")
except KeyError as ke:
    raise KeyError("Please check the column names in the input files.")

try:
    # 输出到一个新的文件，确保输出文件名存在
    if not output_file:
        raise ValueError("Output file path is not provided.")
    filtered_df.to_csv(output_file, index=False)
except ValueError as ve:
    raise ValueError(f"Error while saving the filtered data: {ve}")

if not filtered_df.empty:
    print(f"Filtered data has been saved to {output_file}")
else:
    print("Filtered data is empty; no file was created.")
    raise FileNotFoundError("No data to save to file.")
