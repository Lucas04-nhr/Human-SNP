import os
import pysam
import argparse

def modify_header_and_records(input_vcf, output_vcf):
    with pysam.VariantFile(input_vcf) as vcf_in:
        # Modify the header
        new_header = vcf_in.header.copy()
        for contig in new_header.contigs:
            if contig.startswith('chr'):
                new_contig = contig[3:]  # Remove the 'chr' prefix
                print(f"Modifying header contig from {contig} to {new_contig}")
                new_header.contigs[new_contig] = new_header.contigs[contig]
                del new_header.contigs[contig]


        # Create a VCF Writer with the modified header
        with pysam.VariantFile(output_vcf, 'w', header=new_header) as vcf_out:
            print("Modifying the records...")
            # Modify the records
            for record in vcf_in:
                original_chrom = record.chrom
                if original_chrom.startswith('chr'):
                    new_chrom = original_chrom[3:]  # Remove the 'chr' prefix
                    print(f"Modifying record chromosome from {original_chrom} to {new_chrom}")
                    record.chrom = new_chrom  # Update the record chromosome
            
            # Write the record to the output VCF file
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
print("Input VCF file:\t\t", input_vcf)
print("Output VCF file:\t", output_vcf)
modify_header_and_records(input_vcf, output_vcf)

