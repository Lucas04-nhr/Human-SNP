#!/bin/bash
# Remove unmapped

# Calculate the number of lines
export total_lines=$(wc -l < ./BJ001.head50.sam | grep -oE '[0-9]+')
echo "The total number of lines is $total_lines."
echo "========================================================================"
echo "Processing..."


export i=0

while read -r line; do
    export i=$((i+1))
    echo "Processing line $i..."
    if [[ ! $line =~ \*[\ ]*0[\ ]*0[\ ]*\*[\ ]*\*[\ ]*0[\ ]*0 ]]; then
        echo "$line" >> ./BJ001.head50.removed.sam
        echo "Line $i will be saved."
    else
        echo "Line $i will be removed."
    fi
done < ./BJ001.head50.sam

echo "========================================================================"
echo "Removing unmapped complete."
