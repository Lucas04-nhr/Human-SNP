import pandas as pd
import os
from glob import glob

# Get input and output directories from environment variables
input_dir = os.getenv('INPUT_DIR')
output_dir = os.getenv('OUTPUT_DIR')


# Read cleaned GWAS results file
gwas_df = pd.read_csv(os.path.join(output_dir, "cleaned_gwas_results.csv"))

# Read bacteria allocation file
bacteria_allocation_df = pd.read_csv(os.path.join(output_dir, "bacteria_allocation_full.csv"))
cov_to_species = dict(zip(bacteria_allocation_df["Cov"], bacteria_allocation_df["Species"]))

# Directory containing PLINK result files
plink_results_dir = os.path.join(input_dir)
plink_files = glob(os.path.join(plink_results_dir, "*.assoc.linear"))  # Get all result files

# Result storage dictionary { (SNP, Microbe) -> (BETA, TEST, A1) }
snp_beta_test_a1_dict = {}

# Iterate over all PLINK result files
for i,file in enumerate(plink_files):
    # Parse filename to extract corresponding microbe number
    filename = os.path.basename(file)
    microbe_number = int(filename.split(".")[1][1:])  # Extract number from result.P49.assoc.linear
    microbe_name = cov_to_species.get(microbe_number, None)

    if microbe_name:
        # Read PLINK result file in chunks
        total_lines = sum(1 for line in open(file))
        chunk_size = 10000
        #expected_chunks = (total_lines // chunk_size) + 1
        #print(f"Processing file {i + 1}/{len(plink_files)} with {expected_chunks} chunks")
        for j,chunk in enumerate(pd.read_csv(file, sep='\s+', low_memory=False, chunksize=chunk_size, nrows=500000)):
        # Iterate over PLINK results and store in dictionary
            for _, row in chunk.iterrows():
                snp = row["SNP"]
                beta = row["BETA"]
                test = row["TEST"]
                a1 = row["A1"]
                snp_beta_test_a1_dict[(snp, microbe_name)] = (beta, test, a1)
            # Print progress for each chunk
            print(f"Processed {j + 1} chunks for file {i + 1}/{len(plink_files)}")
    # Print progress
    print(f"Processed {i + 1}/{len(plink_files)} files")

# Create new columns to store BETA, TEST, and A1 values
gwas_df["BETA"] = gwas_df.apply(lambda row: snp_beta_test_a1_dict.get((row["SNP"], row["Cov"]), (None, None, None))[0], axis=1)
gwas_df["TEST"] = gwas_df.apply(lambda row: snp_beta_test_a1_dict.get((row["SNP"], row["Cov"]), (None, None, None))[1], axis=1)
gwas_df["A1"] = gwas_df.apply(lambda row: snp_beta_test_a1_dict.get((row["SNP"], row["Cov"]), (None, None, None))[2], axis=1)

# Save updated GWAS results
gwas_df.to_csv(os.path.join(output_dir,"cleaned_gwas_results_with_beta_test_a1.csv"), index=False)

print("BETA, TEST, and A1 values extracted and saved to cleaned_gwas_results_with_beta_test_a1.csv")