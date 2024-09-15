# 解析 SAM 文件，收集每个 QNAME 对应的 RNEXT 染色体编号。
# 统计每个 QNAME 对应的 RNEXT 染色体出现次数。
# 找到出现次数最多的 RNEXT 染色体编号。
# 仅保留这些 QNAME 的最佳比对行，并输出到新的文件中。

import pysam
from collections import defaultdict, Counter
import argparse

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

def extract_best_rnext(sam_file, output_file):
    # Step 1: Collect RNEXT for each QNAME
    qname_rnext = defaultdict(list)

    with pysam.AlignmentFile(sam_file, "r") as infile:
        for read in infile:
            if read.is_unmapped:  # Skip unmapped reads
                continue
            
            qname = read.query_name  # QNAME
            rnext = read.next_reference_name  # RNEXT

            qname_rnext[qname].append(rnext)

    # Step 2: Find the most common RNEXT for each QNAME
    best_rnext = {}
    for qname, rnext_list in qname_rnext.items():
        counter = Counter(rnext_list)
        most_common_rnext, _ = counter.most_common(1)[0]  # Get the most common RNEXT
        best_rnext[qname] = most_common_rnext

    # Step 3: Write the best RNEXT reads to the output file
    with pysam.AlignmentFile(sam_file, "r") as infile, pysam.AlignmentFile(output_file, "w", header=infile.header) as outfile:
        for read in infile:
            if read.is_unmapped:  # Skip unmapped reads
                continue
            
            qname = read.query_name
            rnext = read.next_reference_name

            # Write the read if the RNEXT is the best RNEXT for the QNAME
            if qname in best_rnext and best_rnext[qname] == rnext:
                outfile.write(read)


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
    exit(127)
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

