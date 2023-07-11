#!/bin/bash
#SBATCH --job-name=SV_LUMPY_TN

module load LUMPY
module load SAMTOOLS

Reference_Genome="$1"
Normal_Bam="$2"
Tumour_Bam="$3"

####################  Normal #####################
# Extract the discordant paired-end alignments.
samtools view -b -F 1294 "$Normal_Bam" > Normal.discordants.unsorted.bam

# Extract the split-read alignments
samtools view -h "$Normal_Bam" \
    | scripts/extractSplitReads_BwaMem -i stdin \
    | samtools view -Sb - \
    > Normal.splitters.unsorted.bam

samtools sort Normal.discordants.unsorted.bam -o Normal.discordants.bam
samtools sort Normal.splitters.unsorted.bam -o Normal.splitters.bam


####################  Tumour #####################
# Extract the discordant paired-end alignments.
samtools view -b -F 1294 "$Tumour_Bam" > Tumour.discordants.unsorted.bam

# Extract the split-read alignments
samtools view -h "$Tumour_Bam" \
    | scripts/extractSplitReads_BwaMem -i stdin \
    | samtools view -Sb - \
    > Tumour.splitters.unsorted.bam

samtools sort Tumour.discordants.unsorted.bam -o Tumour.discordants.bam
samtools sort Tumour.splitters.unsorted.bam -o Tumour.splitters.bam

lumpyexpress \
    -B $Tumour_Bam,$Normal_Bam \
    -S Tumour.splitters.bam,Normal.splitters.bam \
    -D Tumour.discordants.bam,Normal.discordants.bam \
    -o Lumpy.vcf
