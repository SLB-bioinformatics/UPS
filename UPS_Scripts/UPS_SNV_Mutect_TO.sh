#!/bin/bash
#SBATCH --job-name=SNV_Mutect_TO
#SBATCH --mem=10G

module load SAMTOOLS
module load GATK


Reference_Genome="$1"
Tumour_Bam="$2"

if [ ! -f "$Tumour_Bam.bai" ]; then
    samtools index "$Tumour_Bam"
fi

echo "$Reference_Genome"
echo "$Tumour_Bam"

gatk-launch CreateSequenceDictionary -R "$Reference_Genome"

gatk Mutect2  -R "$Reference_Genome" -I "$Tumour_Bam" -O Mutect2_unfiltered.vcf.gz

gatk FilterMutectCalls -R ref.fasta -V Mutect2_unfiltered.vcf.gz -O Mutect2_filtered.vcf