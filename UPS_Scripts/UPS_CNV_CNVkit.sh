#!/bin/bash
#SBATCH --job-name=CNV_CNVkit
#SBATCH --mem=10G
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
        -WE)
            WE="$2"
            shift 2
            ;;
    esac
done



# Print the parsed values for testing
echo "Reference_Genome_Index: $Reference_Genome_Index"
echo "Current_Patient_id: $Current_Patient_id"
echo "TumourFileFlag: $TumourFileFlag"
echo "TumourFile: $TumourFile"
echo "NormalFileFlag: $NormalFileFlag"
echo "NormalFile: $NormalFile"
echo "Support_Path: $Support_Path"
echo "O: $O"
echo "WE: $WE"

if [ "$WE" == "WES" ]; then
    M="hybrid"
    T="$Support_Path/hg38_Twist_ILMN_Exome_2.0_Plus_Panel_annotated.BED"
else
    M="wgs"  
    T=""
fi

refFlat="$Support_Path/refFlat.txt"

# Check for missing required arguments
if [ -z "$Reference_Genome_Index" ] || [ -z "$Current_Patient_id" ]; then
    echo "Error: Missing required arguments."
    exit 1
fi

# Set the paths to the required files
reference_fasta="$1"
tumor_bam="$2"
refFlat="/mnt/beegfs/hassan/Stuart/UniformProcessingOfSamples/TestData/refFlat.txt"

echo "$NormalFile"
echo "$TumourFile"

if [ "$TumourFileFlag" == "Present" ] && [ "$NormalFileFlag" == "Present" ]; then
        cnvkit.py batch "$TumourFile" -n "$NormalFile" -m "$M" -t $T \
                  -f "$reference_fasta" --annotate "$refFlat"

        echo "cnvkit.py batch "$TumourFile" -n "$NormalFile" -m "$M" -t $T \
                  -f "$reference_fasta" --annotate "$refFlat""

elif [ "$TumourFileFlag" == "Present" ] || [ "$NormalFileFlag" == "Present" ]; then
    if [ "$TumourFileFlag" == "Present" ]; then 
        CurrentFile="$TumourFile"
    else
        CurrentFile="$NormalFile"
    fi  
        cnvkit.py batch "$CurrentFile" -n  -m "$M" -t $T \
                  -f "$reference_fasta" --annotate "$refFlat"
        echo "cnvkit.py batch "$CurrentFile" -n  -m "$M" -t $T \
                  -f "$reference_fasta" --annotate "$refFlat""
fi

