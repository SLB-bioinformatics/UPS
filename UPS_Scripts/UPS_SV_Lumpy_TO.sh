#!/bin/bash
#SBATCH --job-name=SV_LUMPY_TO

module load LUMPY
module load SAMTOOLS

Reference_Genome="$1"
Tumour_Bam="$2"


####################  Tumour #####################
# Extract the discordant paired-end alignments.
samtools view -b -F 1294 "$Tumour_Bam" > Tumour.discordants.unsorted.bam

# Extract the split-read alignments
samtools view -h "$Tumour_Bam" \
    | ./lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin \
    | samtools view -Sb - \> Tumour.splitters.unsorted.bam

samtools sort Tumour.discordants.unsorted.bam -o Tumour.discordants.bam
samtools sort Tumour.splitters.unsorted.bam -o Tumour.splitters.bam

lumpyexpress \
    -B $Tumour_Bam \
    -S Tumour.splitters.bam,\
    -D Tumour.discordants.bam, \
    -o Lumpy.vcf
