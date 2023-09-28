#!/bin/bash
# UPS_CNV_CNVkit.sh
#SBATCH --mem=10G
#UPS_CNV_CNVkit_TO.sh
module load CNVKIT

# Initialize variables with default values
Reference_Genome_Index=""
Current_Patient_id=""
Sequencing_Type=""
TumourFileFlag=""
TumourFile=""
NormalFileFlag=""
NormalFile=""


# Loop through command-line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -R)
            Reference_Genome_Index="$2"
            shift 2
            ;;
        -S)
            Current_Patient_id="$2"
            shift 2
            ;;
        -ST)
            Sequencing_Type="$2"
            shift 2
            ;;
        -PTFlaf)
            TumourFileFlag="$2"
            shift 2
            ;;
        -PTFile)
            TumourFile="$2"
            shift 2
            ;;
        -PNFlaf)
            NormalFileFlag="$2"
            shift 2
            ;;
        -PNFile)
            NormalFile="$2"
            shift 2
            ;;
        -Support)
            Support_Path="$2"
            shift 2
            ;;
        -O)
            O="$2"
            shift 2
            ;;
        *)
            echo "Error: Unknown argument $1"
            exit 1
            ;;
    esac
done

# Check for missing required arguments
if [ -z "$Reference_Genome_Index" ] || [ -z "$Current_Patient_id" ] || [ -z "$Sequencing_Type" ]; then
    echo "Error: Missing required arguments."
    exit 1
fi

# Print the parsed values for testing
echo "Reference_Genome_Index: $Reference_Genome_Index"
echo "Current_Patient_id: $Current_Patient_id"
echo "Sequencing_Type: $Sequencing_Type"
echo "TumourFileFlag: $TumourFileFlag"
echo "TumourFile: $TumourFile"
echo "NormalFileFlag: $NormalFileFlag"
echo "NormalFile: $NormalFile"
echo "Support_Path: $Support_Path"
echo "O: $O"

# Set the paths to the required files
reference_fasta="$1"
tumor_bam="$2"
refFlat="/mnt/beegfs/hassan/Stuart/UniformProcessingOfSamples/TestData/refFlat.txt"

echo "$reference_fasta"
echo "$tumor_bam"

if [ "$TumourFileFlag" == "Present" ] && [ "$NormalFileFlag" == "Present" ]; then
    cnvkit.py  "$Tumour_Bam" \
    --normal "$Normal_Bam" \
    --fasta "$Reference_Genome"


else if [ "$TumourFileFlag" == "Present" ] || [ "$NormalFileFlag" == "Present" ]; then
    if [ "$TumourFileFlag" == "Present" ]; then 
        CurrentFile=$TumourFile 
    else
        CurrentFile=$NormalFile 
    fi
    configManta.py \
        cnvkit.py "$CurrentFile" \
        -n  -f "$reference_fasta" \ 
        --annotate "$refFlat"
fi
