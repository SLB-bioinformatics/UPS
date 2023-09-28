#!/bin/bash
#SBATCH --job-name=SNV_Manta_TO
#SBATCH --mem=10G

# Initialize variables with default values
Reference_Genome_Index=""
Current_Patient_id=""
Sequencing_Type=""
TumourFileFlag=""
TumourFile=""
NormalFileFlag=""
NormalFile=""

module load SAMTOOLS
module load GATK


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


module load MANTA

Reference_Genome="$1"
Tumour_Bam="$2"



if [ "$TumourFileFlag" == "Present" ] && [ "$NormalFileFlag" == "Present" ]; then
    ${MANTA_INSTALL_PATH}/bin/configManta.py \
        --normalBam $NormalFile \
        --tumorBam $TumourFile \
        --referenceFasta $Reference_Genome_Index \
        --runDir ${MANTA_ANALYSIS_PATH}

else if [ "$TumourFileFlag" == "Present" ] || [ "$NormalFileFlag" == "Present" ]; then
    if [ "$TumourFileFlag" == "Present" ]; then 
        CurrentFile=$TumourFile 
    else
        CurrentFile=$NormalFile 
    fi
    configManta.py \
        --tumorBam $CurrentFile \
        --referenceFasta $Reference_Genome \
        --runDir .
fi


