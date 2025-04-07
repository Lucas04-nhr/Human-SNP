import os
import pandas as pd
import argparse


# 设置命令行参数解析
parser = argparse.ArgumentParser(description='Filter, add bacteria column, and merge files in a directory.')
parser.add_argument('-d', '--directory', type=str, required=True, help='The directory containing the files to be processed')
parser.add_argument('-o', '--output_file', type=str, required=True, help='The output file path for the combined and filtered results')
#parser.add_argument('-i', '--item', type=str, required=True, help='The item to filter on')
#parser.add_argument('-t', '--threshold', type=float, required=True, help='The threshold value for filtering')
args = parser.parse_args()

# 获取目录路径和输出文件路径
directory = args.directory
output_file = args.output_file
#hreshold = args.threshold

'''
# 设置可排序的列名
sortable_columns = ['P','BETA','A1', 'TEST']

if args.item not in sortable_columns:
    raise ValueError(f"Invalid item to filter on: {args.item}. Please choose from {sortable_columns}")
else:
    sort_columns = args.item
'''

# 检查目录是否存在
if not os.path.isdir(directory):
    raise FileNotFoundError(f"The directory {directory} does not exist.")


# 初始化一个空的DataFrame用于存储合并后的数据
combined_df = pd.DataFrame()

try:
    '''
    # 遍历目录下的所有文件
    for filename in os.listdir(directory):
        if filename.endswith('.assoc.linear'):
            file_path = os.path.join(directory, filename)
            df = pd.read_csv(file_path, sep='\s+')
            print(f"Processing file {filename} with {len(df)} rows and {len(df.columns)} columns.")
            combined_df = pd.concat([combined_df, df], ignore_index=True, sort=False)
    '''
    for filename in os.listdir(directory):
        if filename.endswith('.assoc.linear'):
            file_path = os.path.join(directory, filename)
            if not os.path.isfile(file_path):
                raise FileNotFoundError(f"The file {file_path} does not exist.")

            # 读取文件
            try:
                df = pd.read_csv(file_path, sep='\s+')
            except Exception as e:
                raise ValueError(f"Error reading {file_path}: {e}")

            print(f"Processing file {filename} with {len(df)} rows and {len(df.columns)} columns.")

            # 检查是否包含所需的列
            if 'P' not in df.columns or 'BETA' not in df.columns:
                raise KeyError(f"The file {file_path} does not contain required columns 'P' and 'BETA'.")

            # 筛选出P列和BETA列有值的行
            filtered_df = df.dropna(subset=['P', 'BETA'])
            print(f"Filtered data contains {len(filtered_df)} rows")

            # 获取文件名中的数字并添加细菌列
            try:
                column_number = int(''.join(filter(str.isdigit, filename)))
            except ValueError:
                raise ValueError(f"Error extracting number from filename {filename}.")
            filtered_df.loc[:, 'Cov'] = column_number

            # 合并到总的DataFrame中
            combined_df = pd.concat([combined_df, filtered_df], ignore_index=True, sort=False)

except FileNotFoundError as fnfe:
    raise FileNotFoundError(f"Error while processing files: {fnfe}")
except KeyError as ke:
    raise KeyError(f"Error while processing files: {ke}")
except ValueError as ve:
    raise ValueError(f"Error while processing files: {ve}")

try:
    # 输出到一个新的文件，确保输出文件名存在
    if not output_file:
        raise ValueError("Output file path is not provided.")
    combined_df.to_csv(output_file, index=False)
except ValueError as ve:
    raise ValueError(f"Error while saving the filtered data: {ve}")

if not combined_df.empty:
    print(f"Filtered data has been saved to {output_file}")
else:
    print("Filtered data is empty; no file was created.")
    raise FileNotFoundError("No data to save to file.")