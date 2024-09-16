# 解析 SAM 文件，收集每个 QNAME 对应的 RNEXT 染色体编号。
# 删去插入片段为 0 的行，并输出到新的文件中。

# Test command (execute in the same directory as the script):
# python filter.py -i /mnt/raid6/bacphagenetwork/data/06_unmapped_removed/Beijing/BJ001.removed.sam -o /mnt/raid6/bacphagenetwork/data/08_filtered/Beijing


import pysam
from collections import defaultdict, Counter
import argparse
import sys
import os
import re

# Define exit code
exit_codes = {
    0: "SUCCESS",
    1: "ERROR",
    127: "ERR_NO_SAMPLE_NAME_FOUND"
}

def print_exit_message_and_exit(code):
  if code in exit_codes and code != 0 and code != 1:
    print(f"Exit code: {code} ({exit_codes[code]})")
    sys.exit(code)
  elif code == 0:
    print("The program has run successfully.")
    sys.exit(0)
  elif code == 1:
    print("An uncatched error occurred during the program.")
    sys.exit(1)
  else:
    print("An unknown error occurred during the program.")
    sys.exit(1)

def extract_sample_name(input_file):
    print("Extracting sample name...")
    
    if "GZ" in input_file:
        match = re.search(r'GZ[0-9]{3}', input_file)
    elif "BJ" in input_file:
        match = re.search(r'BJ[0-9]{3}', input_file)

    if match:
        sample_name = match.group(0)
    else:
        sample_name = "NoMatchFound"
    return sample_name

def calculate_lines(sam_file):
    with open(sam_file, "r") as infile:
        lines = 0
        for line in infile:
            lines += 1
    print("Total lines in SAM file: ", lines)
    return lines

def print_log(total_lines, processed_lines):
    percentage = round(processed_lines / total_lines * 100, 2)
    print(f"Processed {processed_lines} lines, ({percentage}%)")

def process_sam_file(infile, best_rnext):

  # Calculate the total number of lines in the SAM file
  total_lines = calculate_lines(infile)
  processed_lines = 0

  # Iterate through the SAM file
  for line in infile:
    # Skip lines with insert size of 0
    if line.template_length == 0:
      continue

    # Get the QNAME and RNEXT
    qname = line.query_name
    rnext = line.next_reference_name
    processed_lines += 1

    # Add the RNEXT to the dictionary
    best_rnext[qname][rnext] += 1

    # Print the progress every 1000000 lines, calculate the percentage of processed lines
    if processed_lines % 1000000 == 0:
      print_log(total_lines, processed_lines)

  # Close the input file
  infile.close()

def write_best_rnext_to_output(input_file, output_file, best_rnext):
  
  # Calculate the total number of lines in the SAM file
  total_lines = calculate_lines(infile)
  processed_lines = 0

  # Iterate through the dictionary to get the best RNEXT for each QNAME
  for qname, rnexts in best_rnext.items():
    best_rnext_value = rnexts.most_common(1)[0][0]

    # Write the best RNEXT to the output file
    infile = pysam.AlignmentFile(input_file, "r")
    for line in infile:
      if line.query_name == qname and line.next_reference_name == best_rnext_value:
        outfile.write(line)
        processed_lines += 1

      # Print the progress every 1000000 lines, calculate the percentage of processed lines
      if processed_lines % 1000000 == 0:
        print_log(total_lines, processed_lines)

    infile.close()

  # Close the output file
  outfile.close()

# MAIN FUNCTION
def extract_best(input_file, output_file):
  # Open the input and output files
  infile = pysam.AlignmentFile(input_file, "r")
  outfile = pysam.AlignmentFile(output_file, "w", template=infile)

  # Initialize a dictionary to store the best RNEXT for each QNAME
  best_rnext = defaultdict(Counter)

  print("Processing SAM file...")
  print("Filtering out insert size of 0...")

  # Process the SAM file
  process_sam_file(infile, best_rnext)

  print("Filtering complete.")
  print("=====================================")
  print("Writing to output file...")

  # Write the best RNEXT to the output file
  write_best_rnext_to_output(input_file, output_file, best_rnext)

  print("Writing complete.")


# Set up the argument parser
parser = argparse.ArgumentParser(description="Process a SAM file to extract the best RNEXT for each QNAME.")
parser.add_argument('-i', '--input', required=True, help="Input SAM file")
parser.add_argument('-o', '--output', required=True, help="Destination of output SAM file")

args = parser.parse_args()

# Get the input and output file names
print("Initializing...")
input_file = args.input
output_path = args.output

# Extract the sample name from the input file
sample_name = extract_sample_name(input_file)

# Check if sample name was extracted
if sample_name == "NoMatchFound":
    print("No sample name found in the input file name. Please check the input file name.")
    print_exit_message_and_exit(127)
else:
    print("Sample name: \t\t" + sample_name)
    output_file = os.path.join(output_path, sample_name + ".best_rnext.sam")

# Print out the input and output file names
print("Input file: \t\t" + input_file)
print("Output file: \t\t" + output_file)
print("Initializing complete.")
print("=====================================")
print("Processing SAM file...")

# Extract the best RNEXT for each QNAME and write to the output file
extract_best(input_file, output_file)

print("Processing complete.")

