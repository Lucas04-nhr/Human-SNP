import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# 设置命令行参数解析
parser = argparse.ArgumentParser(description='Merge and filter files in a directory.')
parser.add_argument('-i', '--input-file', type=str, required=True, help='The file to be processed')
parser.add_argument('-o', '--output-directory', type=str, required=True, help='The output directory for the filtered results')
args = parser.parse_args()

# 获取目录路径和输出文件路径
input_file = args.input_file # /mnt/raid6/bacphagenetwork/data/12_plink_Full/converted/*.assoc.linear
output_directory = args.output_directory # /mnt/raid6/bacphagenetwork/data/12_plink_Full/volcano_plot
if not os.path.exists(output_directory):
  os.makedirs(output_directory)

# 读取表型文件
try:
  df = pd.read_csv(input_file, sep=',', low_memory=False)
  print(f"Input file loaded with {len(df)} rows and {len(df.columns)} columns.")
except FileNotFoundError as fnfe:
  raise FileNotFoundError(f"Error while reading the input file: {fnfe}")

# 构建输出文件名
input_filename = os.path.basename(input_file)
pheno_number = input_filename.split('.')[1][1:].zfill(2)
pheno_name = df['Bacterium'][0]
output_filename = f"{pheno_number}_VolcanoPlot_{pheno_name}.png"
output_file = os.path.join(output_directory, output_filename)
top_csv_file = os.path.join(output_directory, f"00_top_snps.csv")

# 确认是否已存在 top_csv_file
if os.path.exists(top_csv_file):
  print(f"Top SNPs file already exists: {top_csv_file}")
  top_output = pd.read_csv(top_csv_file)
else:
  print(f"Top SNPs file does not exist, creating new one: {top_csv_file}")
  # CHR,SNP,BP,A1,TEST,NMISS,BETA,STAT,P,Bacterium
  top_output = pd.DataFrame(columns=['CHR', 'SNP', 'BP', 'A1', 'TEST', 'NMISS', 'BETA', 'STAT', 'P', 'Bacterium', 'neg_log_p'])

# 数据预处理
# sortable_columns = ['BP', 'NMISS', 'BETA', 'SE', 'R2', 'SIDAK_SD', 'T', 'P']
effective_column = 'BETA'
top_snps = df.nsmallest(int(len(df) * 0.01), 'P')

### TO_DO ###
# 需要考虑此处的top取出来的SNP是否显著
# 暂时就这么输出吧。。。

# 将top_snps添加到top_output中
print(f"Top SNPs selected: {len(top_snps)}")
top_output = pd.concat([top_output, top_snps], ignore_index=True)
top_output.to_csv(top_csv_file, index=False)
print(f"Top SNPs saved to {top_csv_file}")

# # 绘制火山图
# plt.figure(figsize=(10, 6), dpi=100)
# plt.scatter(df[effective_column], df['neg_log_p'], alpha=0.5)
# # 标识显著SNP
# for _, row in top_snps.iterrows():
#   plt.text(row[effective_column], row['neg_log_p'], row['SNP'], fontsize=9)
# plt.axhline(y=-np.log10(0.05), color='r', linestyle='--')
# plt.title(f'Volcano Plot for {pheno_name}')
# plt.xlabel(effective_column)
# plt.ylabel('-log10(P-value)')
# plt.grid()
# plt.savefig(output_file)
# print(f"Volcano plot saved to {output_file}")
