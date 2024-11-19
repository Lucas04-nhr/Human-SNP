import pandas as pd
import os

# 获取要处理的文件列表
input_directory = 'path_to_input_directory'
output_directory = 'path_to_output_directory'
input_files = [f for f in os.listdir(input_directory) if f.endswith('.assoc.linear')]

# 批量处理每个文件
for input_file in input_files:
    # 读取assoc表格文件
    assoc_data = pd.read_csv(os.path.join(input_directory, input_file), delim_whitespace=True)

    # 筛选出test列为add且显著的SNP
    significant_snps = assoc_data[(assoc_data['TEST'] == 'ADD') & (assoc_data['P'] < 0.05)]

    # 保存筛选结果到新的文件，分隔符为\t
    output_file = os.path.join(output_directory, f'significant_{input_file}')
    significant_snps.to_csv(output_file, sep='\t', index=False)