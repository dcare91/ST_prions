#!/bin/bash

# Define paths
orig_fastq_dir="/media/neptun/LocalDisk16TB/Marusa/SpatialTranscr/Raw_Data"
fastq_dir="/home/neptun/src/marusa/spatial_prions/data/fastqs"
meta_file="/home/neptun/src/marusa/spatial_prions/data/metadata/metadata.csv"

# Initialize arrays to store columns
# Initialize arrays for each column
sample_id=()
original_id=()

# Read the CSV file, skipping the header row
header_skipped=false
while IFS=',' read -r col1 col2 _; do
  if [ "$header_skipped" = false ]; then
    header_skipped=true
    continue
  fi
  sample_id+=("$col2")
  original_id+=("$col1")
done < "$meta_file"


# Loop through each pair of sample_id and original_id
for ((i = 0; i < ${#sample_id[@]}; i++)); do
  current_sample_id="${sample_id[i]}"
  current_original_id="${original_id[i]}"
  
  # Construct the find and ln command
  find_command="find \"$orig_fastq_dir/$current_original_id\" -type f -name \"*.fastq.gz\" -exec ln -s {} \"$fastq_dir/$current_sample_id\" \\;"
  
  # Execute the command
  eval "$find_command"
  
  # Check for errors
  if [ $? -ne 0 ]; then
    echo "Error: Failed to create symbolic links for sample_id $current_sample_id and original_id $current_original_id"
  else
    echo "Symbolic links created for sample_id $current_sample_id and original_id $current_original_id"
  fi
done
