#!/bin/bash
#SBATCH --job-name=UPS_Gene_Dirsruption

Sample="$1"

CopyNumberDir="$1/Copy_Number_Variation/CNVkit"
CopyNumberData=$(find "$CopyNumberDir" -type f -name "*.call.cns")
UPS_Gene_Disruption.R "$CopyNumberData"























