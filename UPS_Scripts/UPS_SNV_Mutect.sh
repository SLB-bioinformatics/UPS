#!/bin/bash
#SBATCH --job-name=SNV_Mutect
#SBATCH --mem=10G

module load SAMTOOLS
module load GATK

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
    esac
done

# Check for missing required arguments
if [ -z "$Reference_Genome_Index" ] || [ -z "$Current_Patient_id" ] ; then
    echo "Error: Missing required arguments."
    exit 1
fi

# Print the parsed values for testing
echo "Reference_Genome_Index: $Reference_Genome_Index"
echo "Current_Patient_id: $Current_Patient_id"
echo "TumourFileFlag: $TumourFileFlag"
echo "TumourFile: $TumourFile"
echo "NormalFileFlag: $NormalFileFlag"
echo "NormalFile: $NormalFile"
echo "Support_Path: $Support_Path"
echo "O: $O"



cd "$O/UPS/Samples/$Current_Patient_id/Single_Nucleotide_Variant/Mutect"

mkdir temp


if [ "$TumourFileFlag" == "Present" ] && [ "$NormalFileFlag" == "Present" ]; then
    echo "Performing Tumour-Normal analysis"
    
    gatk Mutect2 \
        --tmp-dir: "temp" \
        -R "$$Reference_Genome_Index" \
        -I "$TumourFile" \
        -I "$NormalFile" \
        -O "Mutect2_unfiltered.vcf.gz"

    gatk GetPileupSummaries \
        -I "$TumourFile" \
        -I "$NormalFile" \
        -V $Support_Path/hg38/small_exac_common_3.hg38.vcf.gz \
        -L $Support_Path/hg38/small_exac_common_3.hg38.vcf.gz \
        -O Getpileupsummaries.table \
        --tmp-dir "temp"

else if [ "$TumourFileFlag" == "Present" ] || [ "$NormalFileFlag" == "Present" ]; then

    if [ "$TumourFileFlag" == "Present" ]; then
        CurSample="$TumourFile"
        echo "Performing Tumour only analysis"
    elif [ "$NormalFileFlag" == "Present" ]; then
        CurSample="$NormalFile"
        echo "Performing Normal only analysis"
    fi

    gatk Mutect2 \
        --tmp-dir:"temp" \
        -R "$$Reference_Genome_Index" \
        -I "$CurSample" \
        -O "Mutect2_unfiltered.vcf.gz"

    gatk GetPileupSummaries \
        -I "$CurSample" \
        -V $Support_Path/hg38/small_exac_common_3.hg38.vcf.gz \
        -L $Support_Path/hg38/small_exac_common_3.hg38.vcf.gz \
        -O Getpileupsummaries.table \
        --tmp-dir "temp"

fi
gatk CalculateContamination \
    -I Getpileupsummaries.table \
    -O Contamination.table \
    --tmp-dir "temp"
    
gatk FilterMutectCalls \
    -R $Reference_Genome \
    -V Tumour_Unfiltered_somatic.vcf.gz \
    --contamination-table Contamination.table \
    -O Filtered_somatic.vcf.gz \
    --tmp-dir "temp"
