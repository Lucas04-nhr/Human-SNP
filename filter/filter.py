# 解析 SAM 文件，收集每个 QNAME 对应的 RNEXT 染色体编号。
# 统计每个 QNAME 对应的 RNEXT 染色体出现次数。
# 找到出现次数最多的 RNEXT 染色体编号。
# 仅保留这些 QNAME 的最佳比对行，并输出到新的文件中。

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

# 1. Collect the RNEXT for each QNAME
def collect_rnext_for_qname(sam_file):
  qname_rnext = defaultdict(list)
  calculate_lines(sam_file)
  with pysam.AlignmentFile(sam_file, "r") as infile:
    for read in infile:
      if read.is_unmapped:  # Jump over the unmapped reads
        continue
      
      qname = read.query_name  # QNAME
      rnext = read.next_reference_name  # RNEXT

      # Log every 100000 reads
      if infile.tell() % 100000 == 0:
          print(f"Processing QNAME: {qname}")
          print(f"RNEXT: \t\t{rnext}")
          # Calculate and print progress percentage
          progress = (infile.tell() / infile.length) * 100
          print(f"Progress: {progress:.2f}%")

      qname_rnext[qname].append(rnext)
  return qname_rnext

# 2. Find the most common RNEXT for each QNAME
def find_most_common_rnext(qname_rnext):
  best_rnext = {}
  for qname, rnext_list in qname_rnext.items():
    counter = Counter(rnext_list)
    most_common_rnext, _ = counter.most_common(1)[0]  # Get the most common RNEXT
    best_rnext[qname] = most_common_rnext
  return best_rnext

# 3. Check if the best RNEXT for each QNAME is all "=" in the RNEXT column
def filter_qname_by_rnext(sam_file, best_rnext):
    qname_validity = defaultdict(bool)

    with pysam.AlignmentFile(sam_file, "r") as infile:
        for read in infile:
            if read.is_unmapped:  # Jump over the unmapped reads
                continue
            
            qname = read.query_name
            rnext = read.next_reference_name

            # If the QNAME is in the best_rnext and the RNEXT is the best_rnext
            if qname in best_rnext and best_rnext[qname] == rnext:
                # Check if the RNEXT is "="
                if rnext == "=":
                    qname_validity[qname] = True  # Set the QNAME to be valid

    return qname_validity

# 4. Write the best RNEXT for each QNAME to the output file
def write_best_rnext_to_output(sam_file, output_file, best_rnext, qname_validity):
  with pysam.AlignmentFile(sam_file, "r") as infile, pysam.AlignmentFile(output_file, "w", header=infile.header) as outfile:
    for read in infile:
      if read.is_unmapped:  # Jump over the unmapped reads
        continue
      
      qname = read.query_name
      rnext = read.next_reference_name

      # If the QNAME is in the best_rnext and the RNEXT is the best_rnext and the QNAME is valid
      if qname in best_rnext and best_rnext[qname] == rnext and qname_validity.get(qname, False):
        outfile.write(read)

def extract_best_rnext(sam_file, output_file):
  # 1. Collect the RNEXT for each QNAME
  qname_rnext = collect_rnext_for_qname(sam_file)

  # 2. Find the most common RNEXT for each QNAME
  best_rnext = find_most_common_rnext(qname_rnext)

  # 3. Check if the best RNEXT for each QNAME is all "=" in the RNEXT column
  qname_validity = filter_qname_by_rnext(sam_file, best_rnext)

  # 4. Write the best RNEXT for each QNAME to the output file
  write_best_rnext_to_output(sam_file, output_file, best_rnext, qname_validity)


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
extract_best_rnext(input_file, output_file)

print("Processing complete.")

