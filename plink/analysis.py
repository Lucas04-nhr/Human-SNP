import gwaslab as gl
import numpy as np
import pandas as pd
import os
import argparse
import re

# Get the input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str, help='Input file path')
parser.add_argument('--output', type=str, help='Output file path')
args = parser.parse_args()

input_path = args.input
output_path = args.output

# Get the list of files
file_list = os.listdir(input_path)
file_list = [f for f in file_list if f.endswith('.assoc.linear')]
chromosome_name = [re.search(r'c(\d{2})', f).group(1) for f in file_list if re.search(r'c(\d{2})', f)]

# Loop through multiple files and draw plots
for file_name in file_list:
  print (f'Processing {file_name}')
  file_path = os.path.join(input_path, file_name)
  
  # Load data
  data = pd.read_table(file_path, sep='\s+')
  
  # Draw the plot
  stats = gl.Sumstats(
    data,
    snpid="SNP",
    chrom="CHR",
    pos="BP",
    p="P",
    verbose=False
  )
  
  # Extract the number from the file name
  i = re.search(r'\d+', file_name).group()
  output_file = f'{output_path}/c{chromosome_name}_{i}.pdf'
  stats.plot_mqq(
    save=output_file,
    save_args={"dpi": 300, "facecolor": "white"}
  )
