import numpy as np
import pandas as pd
import os
import argparse

# Define the function
def adj_merge (input_directory, output_file, nrows_threshold, sort):
  # Get the list of files
  input_files = [os.path.join(input_directory, f) for f in os.listdir(input_directory) if os.path.isfile(os.path.join(input_directory, f))]
  len_files = int(len(input_files)/2)
  # Read the first file
  first_file = input_files[0]
  header = pd.read_csv(first_file, nrows=0, sep='\s+', engine='python')
  # Create an empty DataFrame
  merged_data = pd.DataFrame(columns=header.columns)
  # Read the first five lines of each file and merge
  count = 0
  for input_file in input_files:
    if input_file.endswith(".assoc.linear.adjusted"):
      count += 1
      print('Processing file: ', input_file)
      print('Progress: ', count, '/', len_files)
      data = pd.read_csv(input_file, nrows=nrows_threshold, sep='\s+', engine='python')
      print('Head of the data:')
      print(data.head())
      merged_data = pd.concat([merged_data, data], ignore_index=True)
  # Sort the data
  if sort is not None:
    merged_data = merged_data.sort_values(by=sort)
  # Delete the empty rows
  merged_data = merged_data.dropna(axis=1, how='all')
  # Reset the index
  merged_data.reset_index(drop=True, inplace=True)
  # Move the 3rd row to the last row
  merged_data = pd.concat([merged_data[:2], merged_data[3:], merged_data[2:3]], ignore_index=True)
  # Save the merged data
  merged_data.to_csv(output_file, index=False, sep=',')
  # Clear the memory
  del data
  del merged_data

# Set the arguments
parser = argparse.ArgumentParser()
parser.add_argument('--input-directory', type=str, required=True)
parser.add_argument('--output-file', type=str, required=False)
parser.add_argument('--nrows-threshold', '-n', type=int, required=False, default=5)
parser.add_argument('--sort', '-s', type=str, required=False, default='FDR_BY')

# Parse the arguments
args = parser.parse_args()
input_directory = args.input_directory
output_file = args.output_file
nrows_threshold = args.nrows_threshold
sort = args.sort

# Process the output file name if not given
if output_file is None:
  parent_folder = input_directory.split('/')[-1]
  output_path = os.path.join(input_directory, '..', '..', 'merged')
  # Get the absolute path
  output_path = os.path.abspath(output_path)
  os.makedirs(output_path, exist_ok=True)
  output_file = os.path.join(output_path, f'{parent_folder}.csv')

# Check if the sort parameter is leagal
sort_avail_list = ['CHR', 'SNP', 'UNADJ', 'GC', 'BONF', 'HOLM', 'SIDAK_SS', 'SIDAK_SD', 'FDR_BH', 'FDR_BY', 'Cov']
if sort not in sort_avail_list:
  print('The sort parameter is not available. Please choose from the following list:')
  print(sort_avail_list)
  exit()

# Call the function
print('Merging the files in the directory: ', input_directory)
print("===============================================")
adj_merge(input_directory, output_file, nrows_threshold, sort)
print("===============================================")
print('Output file: ', output_file)
print('Done!')
