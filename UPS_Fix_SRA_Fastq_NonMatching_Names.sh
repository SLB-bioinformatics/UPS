#!/bin/bash

input_file="$1"
Output_file="BamsFixedNames/$1"
# Read the first line of the input file
first_line=$(head -n 1 "$input_file")

# Check if the first line matches the desired formats
if [[ $first_line =~ ^@SRR[0-9]+\.[0-9]+\.[0-9]+$ || $first_line =~ ^@SRR[0-9]+\.[0-9]+/[12]$ ]]; then
  echo "The FASTQ file has read identifiers in the desired formats."
else
    if [ ! -d "BamsFixedNames" ]; then
        mkdir -p "BamsFixedNames"
    fi
    echo "$first_line"
    echo "$Output_file"
  #sbatch check_fastq_format.sh input.fastq
fi
