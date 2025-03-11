import os
import pandas as pd
import argparse
import sys

# 设置命令行参数解析
parser = argparse.ArgumentParser(description='Merge and filter files in a directory.')
parser.add_argument('directory', type=str, help='The directory containing the files to be processed')
parser.add_argument('output_file', type=str, help='The output file path for the filtered results')
args = parser.parse_args()

# 获取目录路径和输出文件路径
directory = args.directory
output_file = args.output_file

# 初始化一个空的DataFrame用于存储合并后的数据
combined_df = pd.DataFrame()

try:
    # 遍历目录下的所有文件
    for filename in os.listdir(directory):
        if filename.endswith('.adjusted'):
            file_path = os.path.join(directory, filename)
            df = pd.read_csv(file_path)
            combined_df = pd.concat([combined_df, df], ignore_index=True)
except Exception as e:
    print(f"Error while processing files: {e}")
    sys.exit(1)

try:
    # 筛选出FDR_BY列小于5*10e-8的行
    filtered_df = combined_df[combined_df['FDR_BY'] <= 5*10e-8]
except Exception as e:
    print(f"Error while filtering data: {e}")
    sys.exit(1)

try:
    # 输出到一个新的文件
    filtered_df.to_csv(output_file, index=False)
except Exception as e:
    print(f"Error while saving the output file: {e}")
    sys.exit(1)

print(f"Filtered data has been saved to {output_file}")
