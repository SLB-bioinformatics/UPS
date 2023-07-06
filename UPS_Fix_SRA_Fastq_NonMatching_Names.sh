#!/bin/bash
module load R
input_file="$1"
UPS_Script_File="$2"
output_file="FastqFixedNames/$(basename "$1")"

# Check if the input file is compressed
if [[ "$input_file" =~ \.gz$ ]]; then
  # Use zcat to decompress the file and read it line by line
  zcat "$input_file" | while IFS= read -r line; do
    # Process each line as needed
    if [[ $line =~ ^@SRR[0-9]+\.[0-9]+\.[0-9]+$ || $line =~ ^@SRR[0-9]+\.[0-9]+/[12]$ ]]; then
      echo "The FASTQ file has read identifiers in the desired formats."
      break
    else
      if [ ! -d "FastqFixedNames" ]; then
        mkdir -p "FastqFixedNames"
      fi
      echo "$line"
      echo "$output_file"
      Rscript "$UPS_Script_File""UPS_Fix_SRA_Fastq_NonMatching_Names.R" "$input_file" "$output_file"
      break
    fi
  done
else
  # Read the first line of the input file
  first_line=$(head -n 1 "$input_file")

  # Check if the first line matches the desired formats
  if [[ $first_line =~ ^@SRR[0-9]+\.[0-9]+\.[0-9]+$ || $first_line =~ ^@SRR[0-9]+\.[0-9]+/[12]$ ]]; then
    echo "The FASTQ file has read identifiers in the desired formats."
  else
    if [ ! -d "FastqFixedNames" ]; then
      mkdir -p "FastqFixedNames"
    fi
    echo "$first_line"
    echo "$output_file"
    Rscript "$UPS_Script_File""UPS_Fix_SRA_Fastq_NonMatching_Names.R" "$input_file" "$output_file"
    fi
fi
