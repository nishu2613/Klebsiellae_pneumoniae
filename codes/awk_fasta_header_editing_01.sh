#!/bin/bash

output_file="all_genomes.faa"

for file in ncbi_dataset/ncbi_dataset/data/GCA*/*protein.faa; do
    gcf_id=$(basename "$(dirname "$file")")
    awk -v gcf="$gcf_id" '/^>/ {print ">" gcf " " substr($0, 2); next} {print}' "$file" >> "$output_file"
done

echo "Merged FASTA saved as: $output_file"