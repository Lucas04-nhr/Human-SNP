import os
import vcf
import argparse

def remove_chr_prefix(input_vcf, output_vcf):
    # Create a VCF Reader
    vcf_reader = vcf.Reader(open(input_vcf, 'r'))
    # Create a VCF Writer
    vcf_writer = vcf.Writer(open(output_vcf, 'w'), vcf_reader)

    for record in vcf_reader:
        print("Processing record:", record)
        print("CHROM field before modification:", record.CHROM)
        # Modify the CHROM field, remove the 'chr' prefix
        if record.CHROM.startswith('chr'):
            print("Removing 'chr' prefix from CHROM field")
            record.CHROM = record.CHROM[3:]  # Remove the 'chr' prefix
            print("CHROM field after modification:", record.CHROM)
        else:
            print("CHROM field does not start with 'chr'")
        
        # Write the record to the output VCF file
        vcf_writer.write_record(record)

    # Close the VCF Reader and Writer
    vcf_writer.close()

# Parse the input arguments
parser = argparse.ArgumentParser(description='Remove the chr prefix from the CHROM field in a VCF file')
parser.add_argument('-i', '--input', required=True, help='Input VCF file')
parser.add_argument('-o', '--output', required=True, help='Output VCF file path')

# Parse the input arguments
args = parser.parse_args()
input_path = args.input
output_path = args.output
input_vcf = os.path.abspath(input_path)
output_vcf = os.path.join(os.path.abspath(output_path), 'modified_' + os.path.basename(input_vcf))

# Call the function to remove the 'chr' prefix
print("Input VCF file:\t", input_vcf)
print("Output VCF file:\t", output_vcf)
remove_chr_prefix(input_vcf, output_vcf)

