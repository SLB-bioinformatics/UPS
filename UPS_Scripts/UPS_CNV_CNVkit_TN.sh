#!/bin/bash
#SBATCH --job-name=CNV_CNVKIT

module load CNVKIT

Reference_Genome="$1"
Normal_Bam="$2"
Tumour_Bam="$3"

cnvkit.py  "$Tumour_Bam" --normal "$Normal_Bam" --fasta "$Reference_Genome"
