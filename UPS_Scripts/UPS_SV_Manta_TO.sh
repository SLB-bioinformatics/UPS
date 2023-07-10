#!/bin/bash
#SBATCH --job-name=SV_MANTA_TO

module load MANTA

Reference_Genome="$1"
Tumour_Bam="$2"

configManta.py \
--tumorBam $Tumour_Bam \
--referenceFasta $Reference_Genome \
--runDir .
