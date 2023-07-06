#!/bin/bash
#SBATCH --job-name=FastqToBam
#SBATCH --mem=10G

module load SAMTOOLS
module load BWA

Reference_Genome="${1}"
File1="${2}"
File2="${3}"
Current_Sample="${4}"

echo "$Reference_Genome"
echo "$File1"
echo "$File2"
echo "$Current_Sample"

# Check if all required arguments are provided
if [ $# -ne 4 ]; then
  echo "Usage: $0 <Reference_Genome> <File1> <File2> <Current_Sample>"
  exit 1
fi

# Aligning the input files
if [ ! -f "${Current_Sample}_1.sai" ]; then
  bwa aln "${Reference_Genome}" "${File1}" > "${Current_Sample}_1.sai"
fi

if [ ! -f "${Current_Sample}_2.sai" ]; then
  bwa aln "${Reference_Genome}" "${File2}" > "${Current_Sample}_2.sai"
fi

# Generating a SAM file by pairing the aligned sequences from both input files with the reference genome
bwa sampe "${Reference_Genome}" "${Current_Sample}_1.sai" "${Current_Sample}_2.sai" "${File1}" "${File2}" > "${Current_Sample}.sam"

# Convert SAM to BAM using Samtools
samtools view -S -b "${Current_Sample}.sam" > "${Current_Sample}.bam"
