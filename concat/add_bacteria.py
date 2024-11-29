import numpy as np
import pandas as pd
import os
import argparse

# Define the function
def add_bacteria_col (input_file, output_file):
  # Read the input file
  data = pd.read_csv(input_file, sep='\t')
  # Get the number of columns
  base_name = os.path.basename(input_file)
  column_number = int(''.join(filter(str.isdigit, base_name)))
  print('Column number: ', column_number)
  # Add the bacteria column
  data['Cov'] = column_number
  # Write the output file
  data.to_csv(output_file, sep='\t', index=False)
  # Clear the memory
  del data


# Set the arguments
parser = argparse.ArgumentParser()
parser.add_argument('--input-file', type=str, required=True)
parser.add_argument('--output-file', type=str, required=False)

# Parse the arguments
args = parser.parse_args()
input_file = args.input_file
output_file = args.output_file

# Process the output file name if not given
if output_file is None:
  input_path = os.path.dirname(input_file)
  parent_folder = input_path.split('/')[-1]
  output_path = os.path.join(input_path, '..', '..', 'modified', parent_folder)
  # Get the absolute path
  output_path = os.path.abspath(output_path)
  os.makedirs(output_path, exist_ok=True)
  output_file = os.path.join(output_path, os.path.basename(input_file))

# Call the function
print('Adding bacteria column to the file: ', input_file)
print('Output file: ', output_file)
add_bacteria_col(input_file, output_file)
print('Done!')
print('')
