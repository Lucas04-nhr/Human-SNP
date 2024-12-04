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
  data['Bacteria'] = column_number
  # Write the output file
  data.to_csv(output_file, sep=',', index=False)
  # Clear the memory
  del data

def set_bacteria_dict (bacteria_file):
  bacteria_dict = pd.Series(bacteria_file.iloc[:, 1].values, index=bacteria_file.iloc[:, 0]).to_dict()
  return bacteria_dict

def replace_bacteria_col (input_file, output_file, bacteria_dict):
  # Read the input file
  data = pd.read_csv(input_file, sep=',')
  # Get the number of columns
  base_name = os.path.basename(input_file)
  column_number = int(''.join(filter(str.isdigit, base_name)))
  print('Column number: ', column_number)
  # Add the bacteria column
  data['Bacteria'] = column_number
  # Replace the bacteria names
  data['Bacteria'] = data['Bacteria'].replace(bacteria_dict)
  # Write the output file
  data.to_csv(output_file, sep=',', index=False)
  # Clear the memory
  del data

# Set the arguments
parser = argparse.ArgumentParser()
parser.add_argument('-a', '--add-bacteria', action='store_true', default=False)
parser.add_argument('-r', '--replace-bacteria', action='store_true', default=False)
parser.add_argument('--input-file', type=str, required=True)
parser.add_argument('--bacteria-file', type=str, required=False)

# Parse the arguments
args = parser.parse_args()
add_bacteria = args.add_bacteria
replace_bacteria = args.replace_bacteria
input_file = args.input_file
bacteria_file = args.bacteria_file

# Process the output file name
input_path = os.path.dirname(input_file)
parent_folder = input_path.split('/')[-1]

# Output path of the modified files
output_path_modified = os.path.join(input_path, '..', '..', 'modified', parent_folder)
# Get the absolute path
output_path_modified = os.path.abspath(output_path_modified)
os.makedirs(output_path_modified, exist_ok=True)
output_file_modified = os.path.join(output_path_modified, os.path.basename(input_file))

# Output path of the replaced files
output_path_replaced = os.path.join(input_path, '..', '..', 'replaced', parent_folder)
# Get the absolute path
output_path_replaced = os.path.abspath(output_path_replaced)
os.makedirs(output_path_replaced, exist_ok=True)
output_file_replaced = os.path.join(output_path_replaced, os.path.basename(input_file))

# Call the function
if add_bacteria:
  print('Adding bacteria column to the file: ', input_file)
  print('Output file: ', output_file_modified)
  add_bacteria_col(input_file, output_file_modified)
  print('Done!')
  print('')
else:
  print('Skipping adding bacteria column.')

if replace_bacteria:
  # Check if the bacteria file is provided
  if not bacteria_file:
    print('Error: bacteria file is not provided.')
    exit()
  # Read the bacteria file
  print('Reading bacteria file: ', bacteria_file)
  bacteria_file = pd.read_csv(bacteria_file)
  # Set the bacteria dictionary
  print('Setting the bacteria dictionary...')
  bacteria_dict = set_bacteria_dict(bacteria_file)
  # Replace the bacteria column
  print('Replacing bacteria column in the file: ', input_file)
  print('Output file: ', output_file_replaced)
  replace_bacteria_col(input_file, output_file_replaced, bacteria_dict)
  print('Done!')
