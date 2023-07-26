#!/bin/bash
#SBATCH --job-name=CNV_CNVKIT_TO
#SBATCH --mem=10G
#UPS_CNV_CNVkit_TO.sh
export MPLCONFIGDIR=$(mktemp -d)

module load CNVKIT
module load BEDTOOLS

Reference_Genome="$1"
Tumour_Bam="$2"
refFlat="/mnt/beegfs/hassan/Stuart/UniformProcessingOfSamples/TestData/refFlat.txt"
echo "$Reference_Genome"
echo "$Tumour_Bam"

# Create a directory to store the CNVkit output
OUTPUT_DIR="$(pwd)"

cnvkit.py batch --method wgs  "$Tumour_Bam" -f "$Reference_Genome" \
    --annotate refFlat.txt --access data/access-5kb-mappable.hg19.bed \
    --output-reference my_flat_reference.cnn

# Perform CNV calling using CNVkit
cnvkit.py batch --method wgs \
  "$Tumour_Bam" \
  --output-dir "$OUTPUT_DIR" \
  --diagram \
  --scatter \
  --reference "$Reference_Genome" \
  --annotate refFlat.txt


