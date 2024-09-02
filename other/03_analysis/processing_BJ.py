import pysam
import argparse
from collections import Counter
import os
import re
import matplotlib.pyplot as plt

def extract_sample_name(input_file):
    match = re.search(r'BJ[0-9]{3}', input_file)
    if match:
        sample_name = match.group(0)
    else:
        sample_name = "NoMatchFound"
    return sample_name

def count_elements_in_tmp_file(tmp_file):
    element_counts = {}

    with open(tmp_file, "r") as infile:
        for line in infile:
            element = line.strip()
            if element in element_counts:
                element_counts[element] += 1
            else:
                element_counts[element] = 1

    counts_list = list(element_counts.values())
    return counts_list

def save_counts_to_static_file(counts_list, static_file):
    with open(static_file, "w") as outfile:
        outfile.write("Flag,Count\n")
        for element, count in element_counts.items():
            outfile.write(f"{element},{count}\n")

def draw_pie_chart(counts_list, output_file):
    plt.figure(figsize=(10, 10), dpi=300)
    plt.pie(counts_list, autopct='%1.1f%%')
    plt.savefig(output_file)


# Set up the argument parser
parser = argparse.ArgumentParser(description="Process a SAM file.")
parser.add_argument('-i', '--input', required=True, help="Input SAM file")
parser.add_argument('-s', '--static', required=False, help="Static file directory")
parser.add_argument('-o', '--output', required=True, help="Output directory")

args = parser.parse_args()

# Get the input and output file names
input_file = args.input
sample_name = extract_sample_name(input_file)
output_dir = args.output
tmp_dir = os.path.join(output_dir, "tmp")
static_dir = args.static
output_file = os.path.join(output_dir, os.path.basename(input_file) + ".piechart.pdf")
tmp_file = os.path.join(tmp_dir, os.path.basename(input_file) + ".tmp")

# Check if sample name was extracted
if sample_name == "NoMatchFound":
    print("No sample name found in the input file name. Please check the input file name.")
    exit(127)

# Create the output directory if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
if not os.path.exists(tmp_dir):
    os.makedirs(tmp_dir)
if not os.path.exists(static_dir):
    os.makedirs(static_dir)

# Open the SAM file
print("Processing file: " + input_file)
samfile = pysam.AlignmentFile(input_file, "r")

# Open the temporary file and write the header
with open(tmp_file, "w") as outfile:
  for read in samfile.fetch():
        outfile.write(f"{read.flag}\n")

# Close the SAM file
samfile.close()

# Open the temporary file and count the flags
counts_list = count_elements_in_tmp_file(tmp_file)

# If -s is provided, save the counts to the static file using *.csv format
if args.static:
    static_file = os.path.join(static_dir, sample_name + ".static.csv")
    save_counts_to_static_file(counts_list, static_file)
    
# Draw the pie chart
draw_pie_chart(counts_list, output_file)

# Remove the temporary file
os.remove(tmp_file)
