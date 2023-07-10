#!/bin/bash
#SBATCH --job-name=SV_LUMPY

module load LUMPY
module load SAMTOOLS

Reference_Genome="$1"
Tumour_Bam="$2"


####################  Tumour #####################
# Extract the discordant paired-end alignments.
samtools view -b -F 1294 "$Tumour_Bam" > Tumour.discordants.unsorted.bam

# Extract the split-read alignments
samtools view -h "$Tumour_Bam" \
    | scripts/extractSplitReads_BwaMem -i stdin \
    | samtools view -Sb - \
    > Tumour.splitters.unsorted.bam

samtools sort Tumour.discordants.unsorted.bam Tumour.discordants.bam
samtools sort Tumour.splitters.unsorted.bam Tumour.splitters.bam

lumpyexpress \
    -B $Tumour_Bam \
    -S Tumour.splitters.bam,\
    -D Tumour.discordants.bam, \
    -o Lumpy.vcf
