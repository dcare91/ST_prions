#!/bin/bash

# define paths and files
proj_dir="/home/neptun/src/marusa/spatial_prions"
meta_file="$proj_dir/data/metadata/metadata.csv"
fastq_main_dir="$proj_dir/data/fastqs"
tiffs_main_dir="$proj_dir/data/images"
jsons_main_dir="$proj_dir/data/manual_alignments"
genome_dir="/home/neptun/src/marusa/Reference_genome/refdata-gex-mm10-2020-A"
out_dir="$proj_dir/output/01st_spaceranger"

cd $out_dir

# Initialize arrays to store the selected columns
original_id=()
sample_id=()
slide_area=()
slide_number=()

# Skip the header
header_skipped=false

# Read the CSV file line by line
while IFS=, read -r -a columns; do
  # Check if this is the header line
  if [ "$header_skipped" = false ]; then
    header_skipped=true
    continue
  fi

  # Extract columns 1, 2, and 5 (0-based index)
  value1="${columns[0]}"
  value2="${columns[1]}"
  value7="${columns[6]}"
  value8="${columns[7]}"

  # Add the values to the respective arrays
  original_id+=("$value1")
  sample_id+=("$value2")
  slide_area+=("$value7")
  slide_number+=("$value8")
done < "$meta_file"

# loop through each element of sample id/original id ...

for ((i = 0; i < ${#sample_id[@]}; i++)); do
  current_sample_id="${sample_id[i]}"
  current_original_id="${original_id[i]}"
  current_slide_area="${slide_area[i]}"
  current_slide_number="${slide_number[i]}"
 
  # look for fastq dir for sample_id
 fastq_dir=$(find $fastq_main_dir -type d -name "$current_sample_id")
 # find manual alignemnt file
 json_f=$(find $jsons_main_dir -type f -name "$current_sample_id*.json")
# find image
 tiff_f=$(find $tiffs_main_dir -type f -name "$current_sample_id*.tif")

echo $current_sample_id
echo $current_original_id
echo $fastq_dir
echo $json_f
echo $tiff_f

  # run spaceranger
  spaceranger count \
    --id sranger_$current_sample_id \
    --transcriptome $genome_dir  \
    --fastqs $fastq_dir \
    --sample $current_original_id \
    --slide $current_slide_number \
    --area $current_slide_area \
    --loupe-alignment $json_f \
    --darkimage $tiff_f \
    --localcores 30 \
    --localmem 57 \
    --no-bam \
    --nosecondary
done
