#!/bin/bash
# Remove duplicates
echo "Removing duplicates from $INPUT_PATH/${sample_name}.marked.sam to $OUTPUT_PATH/${sample_name}.removed.sam..."

while read -r line; do

    if [[ ! $line =~ \*[\ ]*0[\ ]*0[\ ]*\*[\ ]*\*[\ ]*0[\ ]*0 ]]; then
        echo "$line" >> ./BJ001.head50.removed.sam
    fi
done < ./BJ001.head50.sam
