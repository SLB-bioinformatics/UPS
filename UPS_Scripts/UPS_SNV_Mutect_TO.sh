#!/bin/bash
#SBATCH --job-name=SNV_Mutect

module load GATK

Reference_Genome="$1"
Tumour_Bam="$2"

gatk Mutect2 \
     -R "$Reference_Genome" \
     -I "$Tumour_Bam" \
     -O Mutect2.vcf.gz