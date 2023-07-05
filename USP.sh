#!/bin/bash
#SBATCH --job-name=USP
export NCBI_HOME=/mnt/beegfs/hassan/Stuart/
#Extract fastQ from SSR files and prefrom Fastqc on the file to get basic statisics and check sequence quality 

#Set the maximum number of jobs that can run at once
Batch_Job_Limit=5
#Find out how many jobs are currently submitted
Count_Batch_Jobs() {
    Current_Job_Count=$(squeue -u path1357 -h | wc -l)
    #echo=Current_Job_Count
}

#set variables 
Sample_Guide="/mnt/beegfs/hassan/Stuart/UniformProcessingOfSamples/TestData/SampleGuideOneSample.csv"
SSR_Path="/mnt/beegfs/hassan/dbGAP_Ewings/"
Working_Directory_Path="/mnt/beegfs/hassan/Stuart/UniformProcessingOfSamples/"
Reference_Genome="/mnt/beegfs/hassan/Stuart/UniformProcessingOfSamples/TestData/Index/Homo_sapiens_assembly38.fasta"
UPS_Script_File="/mnt/beegfs/hassan/Stuart/UniformProcessingOfSamples/Scripts/UPS_Scripts/"
Chromosome_Of_Intrest="Chr17"
Gene_of_Intrest="TP53"


Count_Batch_Jobs() {
squeue -u path1357 -t PD -o %i | wc -l
}

# Read the CSV file and extract the specified column
Patient_ids=$(cut -d ',' -f 4 "$Sample_Guide")



# Create a list of unique values from the column
Unique_Patient_ids=$(echo "$Patient_ids" | awk -F, 'NR>1 && !seen[$1]++ {print $1}')

#create USP dir
if [ ! -d "$Working_Directory_Path""UPS" ]; then
    mkdir -p "$Working_Directory_Path""UPS"
fi

#create USP dir
if [ ! -d "$Working_Directory_Path""UPS/Samples" ]; then
    mkdir -p "$Working_Directory_Path""UPS/Samples"
fi

#Create Sample folders and Exstract Fastq from sra files
for Current_Patient_id in $Unique_Patient_ids; do
    #Create a folder for each Patient 
    if [ ! -d "$Working_Directory_Path""UPS/Samples/""$Current_Patient_id" ]; then
        mkdir -p "$Working_Directory_Path""UPS/Samples/""$Current_Patient_id"
    fi

    #Create a folder for each Patient FastQ
    if [ ! -d "$Working_Directory_Path""UPS/Samples/""$Current_Patient_id/FastQ" ]; then
        mkdir -p "$Working_Directory_Path""UPS/Samples/$Current_Patient_id/FastQ"
    fi

    echo "Processing Patient ID: $Current_Patient_id"
    # Find rows that match the current value in column 4
    matching_rows=$(awk -F ',' -v col4="$Current_Patient_id" 'BEGIN {OFS=","} $4 == col4' "$Sample_Guide")
    # Select the matching rows from column 5
    matching_values=$(echo "$matching_rows" | cut -d ',' -f 8| tr -d '\r')
    for Current_SSR in $matching_values; do
        echo "$SSR_Path""$Current_SSR"
        cd "$Working_Directory_Path/UPS/Samples/$Current_Patient_id/FastQ"
        if [ ! -f "$Current_SSR""_pass_1.fastq.gz" ]; then
            #Extract FastQ from SSR file
            echo "Extraction Needed for $Current_SSR"
            while [[ $(Count_Batch_Jobs) -gt $Batch_Job_Limit ]]
            do
                sleep 60
            done
            File_For_Extraction="$SSR_Path""$Current_SSR"
            echo "$File_For_Extraction"
            sbatch "$UPS_Script_File""UPS_Fastq-dump.sh" "$File_For_Extraction"

        fi
        # check file names are in the correct format
        sbatch "$UPS_Fix_SRA_Fastq_NonMatching_Names.sh" "${Current_SSR}_pass_1.fastq.gz"
    done
done 

#Wait untill all Fastq_dump jobs have finished running 
Fastq_Dump_Jobs_Running=$(squeue -u path1357 -h -o "%i" | grep "Fastq-dump" | wc -l)
while [ $Fastq_Dump_Jobs_Running -gt 0 ]
do
    sleep 60
    Fastq_Dump_Jobs_Running=$(squeue -u path1357 -h -o "%i" | grep "Fastq-dump" | wc -l)
done

#Create Bams 
for Current_Patient_id in $Unique_Patient_ids; do
    echo "Processing Patient ID: $Current_Patient_id"
    matching_rows=$(awk -F ',' -v col4="$Current_Patient_id" 'BEGIN {OFS=","} $4 == col4' "$Sample_Guide")
    matching_values=$(echo "$matching_rows" | cut -d ',' -f 8 | tr -d '\r')
    
    #Create a folder for each Patient bam
    if [ ! -d "$Working_Directory_Path""UPS/Samples/""$Current_Patient_id/BAM" ]; then
        mkdir -p "$Working_Directory_Path""UPS/Samples/$Current_Patient_id/BAM"
    fi
    cd "$Working_Directory_Path""UPS/Samples/$Current_Patient_id/BAM"
    
    #Create convert Fastq to Bam by mapping to reference genome
    for Current_SSR in $matching_values; do
        if [ ! -f "${Current_SSR}.bam" ]; then
            echo "BAM Needed for $Current_SSR"
            while [[ $(Count_Batch_Jobs) -gt $Batch_Job_Limit ]]; do
                sleep 60
            done
            File1="${Working_Directory_Path}UPS/Samples/$Current_Patient_id/FastQ/${Current_SSR}_pass_1.fastq.gz"
            File2="${Working_Directory_Path}UPS/Samples/$Current_Patient_id/FastQ/${Current_SSR}_pass_2.fastq.gz"
            #sbatch "$UPS_Script_File""UPS_Fastq_To_Bam.sh" "$Reference_Genome" "$File1" "$File2" "$Current_SSR"
        fi
    done
done


#Wait untill all Fastq_To_Bam have finished running Hashed out
#Fastq_Dump_Jobs_Running=$(squeue -u path1357 -h -o "%i" | grep "Fastq-dump" | wc -l)
#while [ $Fastq_Dump_Jobs_Running -gt 0 ]
#do
#    sleep 60
#    Fastq_Dump_Jobs_Running=$(squeue -u path1357 -h -o "%i" | grep "Fastq-dump" | wc -l)
#done

# calls SNVs SV and CNVs
for Current_Patient_id in $Unique_Patient_ids; do
    #Create a folder for each Patient SNV
    if [ ! -d "$Working_Directory_Path""UPS/Samples/""$Current_Patient_id/Single_Nucleotide_Variant" ]; then
        mkdir -p "$Working_Directory_Path""UPS/Samples/$Current_Patient_id/Single_Nucleotide_Variant"
    fi

    #Create a folder for each Patient SV
    if [ ! -d "$Working_Directory_Path""UPS/Samples/""$Current_Patient_id/Structural_Variant" ]; then
        mkdir -p "$Working_Directory_Path""UPS/Samples/$Current_Patient_id/Structural_Variant"
    fi

    #Create a folder for each Patient CNV
    if [ ! -d "$Working_Directory_Path""UPS/Samples/""$Current_Patient_id/Copy_Number_Variation" ]; then
        mkdir -p "$Working_Directory_Path""UPS/Samples/$Current_Patient_id/Copy_Number_Variation"
    fi

    echo "Processing Patient ID: $Current_Patient_id"
    Patient_Tumour_WGS_SSR=$(awk -F ',' -v col4="$Current_Patient_id" 'BEGIN {OFS=","} $4 == col4 && $6 == "Tumor" && $7 == "WGS" {print $8}' "$Sample_Guide" | tr -d '\r')
    Patient_Normal_WGS_SSR=$(awk -F ',' -v col4="$Current_Patient_id" 'BEGIN {OFS=","} $4 == col4 && $6 == "Normal" && $7 == "WGS" {print $8}' "$Sample_Guide" | tr -d '\r')
    Patient_Tumour_WGS_Bam="$Working_Directory_Path""UPS/Samples/$Current_Patient_id/BAM/$Patient_Tumour_WGS_SSR.bam"
    Patient_Normal_WGS_Bam="$Working_Directory_Path""UPS/Samples/$Current_Patient_id/BAM/$Patient_Normal_WGS_SSR.bam"
    echo "$Patient_Tumour_WGS_Bam"
    echo "$Patient_Normal_WGS_Bam"

    if [[ -n "$Patient_Tumour_WGS_SSR" && -n "$Patient_Normal_WGS_SSR" ]]; then
        echo "$Current_Patient_id Tumour vs Normal analysis "
    
        ##########################   Calling SNVs    ##########################
        #Create a folder for Mutec
        if [ ! -d "$Working_Directory_Path""UPS/Samples/""$Current_Patient_id/Single_Nucleotide_Variant/Mutect" ]; then
            mkdir -p "$Working_Directory_Path""UPS/Samples/$Current_Patient_id/Single_Nucleotide_Variant/Mutect"
        fi

        #run Mutect
        cd "$Working_Directory_Path""UPS/Samples/$Current_Patient_id/Single_Nucleotide_Variant/Mutect"
        if [ ! -f "tumorSV.vcf.gz" ]; then
            while [[ $(Count_Batch_Jobs) -gt $Batch_Job_Limit ]]; do
                    sleep 60
            done
            #sbatch "$UPS_Script_File""UPS_SV_Mutect.sh" "$Reference_Genome" "$Patient_Normal_WGS_Bam" "$Patient_Tumour_WGS_Bam" 
        fi

        ##########################   Calling SVs    ##########################

        #run manta
        if [ ! -d "$Working_Directory_Path""UPS/Samples/""$Current_Patient_id/Structural_Variant/Manta" ]; then
            mkdir -p "$Working_Directory_Path""UPS/Samples/$Current_Patient_id/Structural_Variant/Manta"
        fi

        cd "$Working_Directory_Path""UPS/Samples/$Current_Patient_id/Structural_Variant/Manta"
        if [ ! -f "tumorSV.vcf.gz" ]; then
            while [[ $(Count_Batch_Jobs) -gt $Batch_Job_Limit ]]; do
                    sleep 60
            done
            #sbatch "$UPS_Script_File""UPS_SV_Manta.sh" "$Reference_Genome" "$Patient_Normal_WGS_Bam" "$Patient_Tumour_WGS_Bam" 
        fi

        #run Lumpy
        if [ ! -d "$Working_Directory_Path""UPS/Samples/""$Current_Patient_id/Structural_Variant/Lumpy" ]; then
            mkdir -p "$Working_Directory_Path""UPS/Samples/$Current_Patient_id/Structural_Variant/Lumpy"
        fi

        cd "$Working_Directory_Path""UPS/Samples/$Current_Patient_id/Structural_Variant/Lumpy"
        if [ ! -f "Lumpy.vcf" ]; then
            while [[ $(Count_Batch_Jobs) -gt $Batch_Job_Limit ]]; do
                    sleep 60
            done
            #sbatch "$UPS_Script_File""UPS_SV_Lumpy.sh" "$Reference_Genome" "$Patient_Normal_WGS_Bam" "$Patient_Tumour_WGS_Bam" 
        fi

        ##########################   Calling CNVs    ##########################
        #Create a folder for Manta
        if [ ! -d "$Working_Directory_Path""UPS/Samples/""$Current_Patient_id/Copy_Number_Variation/CNVkit" ]; then
            mkdir -p "$Working_Directory_Path""UPS/Samples/$Current_Patient_id/Copy_Number_Variation/CNVkit"
        fi

        #run manta
        cd "$Working_Directory_Path""UPS/Samples/$Current_Patient_id/Copy_Number_Variation/CNVkit"
        if [ ! -f "tumorSV.vcf.gz" ]; then
            while [[ $(Count_Batch_Jobs) -gt $Batch_Job_Limit ]]; do
                    sleep 60
            done
            #sbatch "$UPS_Script_File""UPS_CNV_CNVKit.sh" "$Reference_Genome" "$Patient_Normal_WGS_Bam" "$Patient_Tumour_WGS_Bam" 
        fi
    elif [[ -n "$Patient_Tumour_WGS_SSR" ]]; then
        echo "$Current_Patient_id Tumour only analysis "
    else
        echo "Error: No genomic information found for $Current_Patient_id"
    fi
done 