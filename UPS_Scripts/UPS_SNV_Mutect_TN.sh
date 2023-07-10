#!/bin/bash
#SBATCH --job-name=SNV_Mutect_TN

module load GATK

Reference_Genome="$1"
Normal_Bam="$2"
Tumour_Bam="$3"

gatk Mutect2 \
     -R "$Reference_Genome" \
     -I "$Tumour_Bam" \
     -I "$Normal_Bam" \
     -O Mutect2.vcf.gz