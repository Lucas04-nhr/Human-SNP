#! /bin/bash

# Change to the directory where the files are located
cd /mnt/raid6/bacphagenetwork/data/bwa_analysis

# List all files in the folder ending with .sam
for file in GZ*.sam; do
    # Check if the file name matches the pattern
    if ! [[ $file =~ ^GZ[0-9]{3}\.sam$ ]]; then
        # Extract the number from the file name
        number=$(echo $file | sed -E 's/GZ([0-9]+)\.sam/\1/')
        # Create the new name
        new_name=$(printf "GZ%03d.sam" $number)
        # Rename the file
        mv "$file" "$new_name"
        echo "Renamed $file to $new_name"
    fi
done
