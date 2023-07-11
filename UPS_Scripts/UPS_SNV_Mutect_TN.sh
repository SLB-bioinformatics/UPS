#!/bin/bash
#SBATCH --job-name=SNV_Mutect_TN
#SBATCH --mem=10G

module load GATK
Reference_Genome="$1"
Normal_Bam="$2"
Tumour_Bam="$3"

if [ ! -f "$Normal_Bam.bai" ]; then
    samtools index "$Normal_Bam"
fi

if [ ! -f "$Tumour_Bam.bai" ]; then
    samtools index "$Tumour_Bam"
fi

gatk Mutect2  -R "$Reference_Genome"  -I "$Tumour_Bam" -I "$Normal_Bam" -O Mutect2.vcf.gz