#!/bin/bash
#SBATCH --job-name=CNV_CNVKIT_TO

export MPLCONFIGDIR=$(mktemp -d)

module load CNVKIT

Reference_Genome="$1"
Tumour_Bam="$2"

cnvkit.py  "$Tumour_Bam" --normal "$Normal_Bam" --fasta "$Reference_Genome"
