#!/bin/bash
#SBATCH --job-name=SV_MANTA

module load MANTA

Reference_Genome="$1"
Normal_Bam="$2"
Tumour_Bam="$3"

configManta.py \
--normalBam $Normal_Bam \
--tumorBam $Tumour_Bam \
--referenceFasta $Reference_Genome \
--runDir .