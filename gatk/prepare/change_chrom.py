import os
import pysam
import argparse

def modify_chrom_prefix(input_vcf, output_vcf):
    with pysam.VariantFile(input_vcf) as vcf_in, pysam.VariantFile(output_vcf, 'w', header=vcf_in.header) as vcf_out:
        for record in vcf_in:
            print("CHROM field before modification:", record.chrom)
            # Modify the CHROM field, remove the 'chr' prefix
            if record.chrom.startswith('chr'):
                print("Removing 'chr' prefix from CHROM field")
                record.chrom = record.chrom[3:]  # Remove the 'chr' prefix
                print("CHROM field after modification:", record.chrom)
            else:
                print("CHROM field does not start with 'chr'")
            
            # Write the modified record to the output VCF file
            vcf_out.write(record)


# Parse the input arguments
parser = argparse.ArgumentParser(description='Remove the chr prefix from the CHROM field in a VCF file')
parser.add_argument('-i', '--input', required=True, help='Input VCF file')
parser.add_argument('-o', '--output', required=True, help='Output VCF file path')

# python change_chrom.py -i /mnt/raid6/bacphagenetwork/data/00_bwa_index/GRCh38/known-sites/Homo_sapiens_assembly38.dbsnp138.vcf -o /mnt/raid6/bacphagenetwork/data/00_bwa_index/GRCh38/known-sites/

# Parse the input arguments
args = parser.parse_args()
input_path = args.input
output_path = args.output
input_vcf = os.path.abspath(input_path)
output_vcf = os.path.join(os.path.abspath(output_path), 'modified_' + os.path.basename(input_vcf))

# Call the function to remove the 'chr' prefix
print("Input VCF file:\t", input_vcf)
print("Output VCF file:\t", output_vcf)
modify_chrom_prefix(input_vcf, output_vcf)

