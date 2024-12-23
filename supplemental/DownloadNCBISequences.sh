#!/bin/bash

# File containing the list of IDs
id_file="Lu2018_BatKobuvirusSequences.txt"

# Output directory for the FASTA files
output_dir="fasta_files"
mkdir -p "$output_dir"

# Loop over each ID in the file and run the efetch command
while IFS= read -r id; do
    echo "Fetching sequence for ID: $id"
    efetch -db nucleotide -id "$id" -format fasta > "$output_dir/$id.fasta"
done < "$id_file"

echo "All sequences have been fetched and saved to $output_dir"