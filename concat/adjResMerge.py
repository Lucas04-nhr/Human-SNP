import os
import pandas as pd

# 设置文件路径
input_directory = '/mnt/raid6/bacphagenetwork/data/12_plink/Beijing/results/modified/c004'
output_file = '/mnt/raid6/bacphagenetwork/data/12_plink/Beijing/results/modified/c004/merged_first_five_lines.csv'
input_files = [os.path.join(input_directory, f) for f in os.listdir(input_directory) if os.path.isfile(os.path.join(input_directory, f))]

# 读取第一个文件的表头
first_file = input_files[0]
header = pd.read_csv(first_file, nrows=0)

# 初始化一个空的 DataFrame 来存储合并的结果，并添加表头
merged_data = pd.DataFrame(columns=header.columns)

# 读取每个文件的前五行（不包含表头）并合并
for input_file in input_files:
    data = pd.read_csv(input_file, skiprows=1, nrows=5)
    merged_data = pd.concat([merged_data, data], ignore_index=True)

# 保存合并结果到新的文件
merged_data.to_csv(output_file, index=False)