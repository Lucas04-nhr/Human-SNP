import gwaslab as gl
import numpy as np
import pandas as pd
import os
import argparse
import re

# Tell the user what is happening
print("THIS IS A TEST SCRIPT")
# Wait for 3 seconds untill the script continues
time.sleep(3)

# Get the input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str, help='Input file')
parser.add_argument('--output', type=str, help='Output file path')
args = parser.parse_args()

input_file = args.input
output_path = args.output

# Get the list of files
input_path = os.path.dirname(input_file)
file_name = input_file

print (f'Processing {file_name}')
file_path = os.path.join(input_path, file_name)

# Load data
stats = pd.read_table(file_path, sep='\s+', low_memory=False)

# Sort the data
stats.sort_values(by="P", inplace=True)

# Set the threshold p-value
threshold_pvalue = 5e-8

# Get the chromosomes by unique values
chromosomes = stats['CHR'].unique()
highlight_snps = []

for chrom in chromosomes:
    chrom_data = stats[stats['CHR'] == chrom]
    significant_snps = chrom_data[chrom_data['P'] < threshold_pvalue]
    total_snps = len(chrom_data)

    # Calculate the percentage of significant SNPs
    if total_snps > 0:
        percentage = len(significant_snps) / total_snps
    else:
        percentage = 0

    # If the percentage is greater than 0.5, highlight the most significant SNP
    if percentage >= 0.5:
        if not significant_snps.empty:
            most_significant_snp = significant_snps.iloc[0]
            highlight_snps.append(most_significant_snp['SNP'])

# Draw the plot
output_file = f'{output_path}/c005_test.pdf'
stats.plot_mqq(
    highlight=highlight_snps,  # Specify the SNPs to highlight
    highlight_windowkb=500,  # Specify the window size for highlighting
    highlight_color="#CB132D",  # Specify the color for highlighting
    anno=True,  # Add annotations
    anno_fontsize=8,  # Set the font size for annotations
    save=output_file,
    font_family="DejaVu Sans",
    save_args={"dpi": 300, "facecolor": "white"}
)

# Draw the plot
# stats = gl.Sumstats(
#   data,
#   snpid="SNP",
#   chrom="CHR",
#   pos="BP",
#   p="P",
#   verbose=False
# )

# Extract the number from the file name
# i = re.search(r'\d+', file_name).group()
# output_file = f'{output_path}/c005_test.pdf'
# stats.plot_mqq(
#   save=output_file,
#   font_family="DejaVu Sans",
#   save_args={"dpi": 300, "facecolor": "white"}
# )

# Clean up
del stats
