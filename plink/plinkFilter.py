import os
import pandas as pd

# 指定目录路径
directory = 'data/12_plink/Full/modified/bac_age'

# 初始化一个空的DataFrame用于存储合并后的数据
combined_df = pd.DataFrame()

# 遍历目录下的所有文件
for filename in os.listdir(directory):
    if filename.endswith('.adjusted'):
        file_path = os.path.join(directory, filename)
        df = pd.read_csv(file_path)
        combined_df = pd.concat([combined_df, df], ignore_index=True)

# 筛选出FDR_BY列小于5*10e-8的行
filtered_df = combined_df[combined_df['FDR_BY'] <= 5*10e-8]

# 输出到一个新的文件
output_file = 'Full/merged/filtered_bac_age.csv'
filtered_df.to_csv(output_file, index=False)

print(f"Filtered data has been saved to {output_file}")