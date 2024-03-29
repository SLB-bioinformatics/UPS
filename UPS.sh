#!/bin/bash
#SBATCH --job-name=USP_Orchestrator

# sbatch UPS/UPS.sh -SG SampleMapOSSupperShort.csv -CF Configuration_File.YAML  -O Output/

i=$(pwd)

export NCBI_HOME=/mnt/beegfs/hassan/Stuart/


####################################
# Functions
####################################

# Find out how many total jobs are currently submitted by the user.
Count_Total_Jobs() {
    Total_Job_Count=$(squeue -u $(whoami) -h | wc -l)
}

# Find out how many total jobs submitted by the user are currently running.
Count_Job_Type() {
    local Job_Name="$1"
    Current_Job_Count=$(squeue -u $(whoami) -h | grep "$Job_Name" | wc -l)
}

#################################### 
# Parses Command Line 
#################################### 

sg_flag=false  
cf_flag=false  
O="." # Sets current working directory to output directory by defualt
for arg in "$@"; do # Parses the command line input to identife Sample Map, Configuration YAML and output directory
    case "$arg" in
        -SG)
            if [[ $# -gt 1 && ! "$2" =~ ^- ]]; then 
                sg_argument="$2"  
                sg_flag=true 
                shift 2 
            else
                echo "Error: -SG argument requires a a Sample guide .csv"
                exit 1
            fi
            ;;
        -CF)
            if  [[ $# -gt 1 && ! "$2" =~ ^- ]]; then
                cf_argument="$2" 
                cf_flag=true  
                shift 2 
            else
                echo "Error: -CF argument requires a configuation file"
                exit 1
            fi
            ;;
        -O)
            if  [[ $# -gt 1 && ! "$2" =~ ^- ]]; then
                O="$2"  
                shift 2 
            fi
            ;;
    esac
done

#################################### 
# Get Setting From Configuration File YAML 
####################################

while IFS= read -r line; do
    # Use case statements to extract values based on keys
    case "$line" in
        Batch_Job_Limit:*) Batch_Job_Limit="${line#*: }";;
        Reference_Genome:*) Reference_Genome="${line#*: }";;
        UPS_Script_File:*) UPS_Script_File="${line#*: }";;
        Gene_of_Intrest:*) Gene_of_Intrest="${line#*: }";;
        Packages:*) Packages="${line#*: }";;
        Support_Path:*) Support_Path="${line#*: }";;
    esac
done < "$cf_argument"

if    [ -z "$Batch_Job_Limit" ] || [ -z "$Reference_Genome" ] || [ -z "$UPS_Script_File" ] || [ -z "$Packages" ] ||  [ -z "$Support_Path" ]; then
  echo "Error! Configuration File YAML validation failed missing one or more inputs.
        Requires
        Batch_Job_Limit:    The max number of jobs to run at once.
        Reference_Genome:   Path to the reference genome fasta file.
        UPS_Script_File:    Path to where the UPS scripts are located.
        Gene_of_Intrest:    Not required but priotisies plots for this gene in the report.
        Packages:           Path to where the required Packages are located.
        Support_Path:       Path to where the UPS support data"
  exit 1
fi


cd $O
O=$(pwd)


#################################### 
# Validate Sample Map csv
####################################

echo "Validating Sample Map"
## hashed out for development on non slurm systems 
#Rscript "${UPS_Script_File}UPS_Sample_Map_Validation.R" $sg_argument

#while ! [ -e "Sample_Map_Validation_Report.csv" ]; do
#  sleep 10
#done


#N_invalid=$(awk -F',' 'inlsdex($4, "Invalid") {count++} END {print count}' Sample_Map_Validation_Report.csv)
#if [[ "$N_invalid" -gt 0 ]]; then
#    echo "Error! Sample map validation failed missing coloumns or incorrect data. 
#Please see Sample_Map_Validation_Report.csv for more details. "
#    exit 1
#fi

#################################### 
# Set Up Directory Structure 
####################################

#SM_array="Sample_Map_Validation_Report.csv"
#Sample_ID_pos=$(awk -F',' 'NR==2 { print $3 }' "$SM_array")
#Cohort_pos=$(awk -F',' 'NR==3 { print $3 }' "$SM_array")
#Patient_ID_pos=$(awk -F',' 'NR==4 { print $3 }' "$SM_array")
#Cancer_Type_pos=$(awk -F',' 'NR==5 { print $3 }' "$SM_array")
#Sequencing_Type_pos=$(awk -F',' 'NR==6 { print $3 }' "$SM_array")
#Data_Format_pos=$(awk -F',' 'NR==7 { print $3 }' "$SM_array")
#Source_pos=$(awk -F',' 'NR==8 { print $3 }' "$SM_array")
#Sex_pos=$(awk -F',' 'NR==9 { print $3 }' "$SM_array")
#Location_Of_Data_pos=$(awk -F',' 'NR==10 { print $3 }' "$SM_array")

#create USP and directorys 
if [ ! -d "$O""/UPS" ]; then
    mkdir -p "$O""/UPS"
fi

#create USP dir
if [ ! -d "$O""/UPS/Samples" ]; then
    mkdir -p "$O""/UPS/Samples"
fi

if [ ! -d "$O""/UPS/Reports" ]; then
    mkdir -p "$O""/UPS/Reports"
fi

# Read the CSV file line by line and split each line into fields
while IFS=',' read -r Sample_ID Cohort Patient_ID Cancer_Type Sequencing_Type Data_Format Source Sex Location_Of_Data
do  
    SM_array+=("$Sample_ID, $Cohort, $Patient_ID, $Cancer_Type, $Sequencing_Type, $Data_Format, $Source, $Sex, $Location_Of_Data")
done < "$i/$sg_argument"

for element in "${SM_array[@]}"
do
    # Use awk to extract the Patient_ID field (assuming it's the third field)
    patient_id=$(echo "$element" | awk -F', ' '{print $3}')
    # Check if the patient_id is not empty and is not already in the unique_patient_ids array
    if [ -n "$patient_id" ] && [ "$patient_id" != "Patient_ID" ] && ! [[ " ${unique_patient_ids[@]} " =~ "$patient_id" ]]; then
        # Add the unique patient_id to the array
        unique_patient_ids+=("$patient_id")
        
    fi
done

#set up folder system for samples
for Current_Patient_id in "${unique_patient_ids[@]}"
do
    #Create a folder for each Patient 
    if [ ! -d "$O""/UPS/Samples/""$Current_Patient_id" ]; then
        mkdir -p "$O""/UPS/Samples/""$Current_Patient_id"
    fi

    if [ ! -d "$O""/UPS/Reports/""$Current_Patient_id" ]; then
        mkdir -p "$O""/UPS/Reports/""$Current_Patient_id"
    fi

    for line in "${SM_array[@]}"; do
        if [[ "$line" == *"$Current_Patient_id"*"SRR"* ]]; then
            if [ ! -d "$O/UPS/Samples/$Current_Patient_id/FastQ" ]; then
                mkdir -p "$O/UPS/Samples/$Current_Patient_id/FastQ"
            fi
            break
        fi

        if [[ "$line" == *"$Current_Patient_id"*"SRR"* ]] || [[ "$line" == *"$Current_Patient_id"*"FastQ"* ]]; then
            if [ ! -d "$O/UPS/Samples/$Current_Patient_id/FastQ" ]; then
                mkdir -p "$O/UPS/Samples/$Current_Patient_id/FastQ"
            fi
            break
        fi
    done

    if [ ! -d "$O/UPS/Samples/""$Current_Patient_id/BAM" ]; then
        mkdir -p "$O/UPS/Samples/$Current_Patient_id/BAM"
    fi

    if [ ! -d "$O/UPS/Samples/""$Current_Patient_id/Single_Nucleotide_Variant" ]; then
        mkdir -p "$O/UPS/Samples/$Current_Patient_id/Single_Nucleotide_Variant"
    fi

     if [ ! -d "$O/UPS/Samples/""$Current_Patient_id/Single_Nucleotide_Variant/Mutect" ]; then
        mkdir -p "$O/UPS/Samples/$Current_Patient_id/Single_Nucleotide_Variant/Mutect"
    fi

    if [ ! -d "$O/UPS/Samples/""$Current_Patient_id/Structural_Variant" ]; then
        mkdir -p "$O/UPS/Samples/$Current_Patient_id/Structural_Variant"
    fi

    if [ ! -d "$O/UPS/Samples/""$Current_Patient_id/Structural_Variant/Manta" ]; then
        mkdir -p "$O/UPS/Samples/$Current_Patient_id/Structural_Variant/Manta"
    fi

    if [ ! -d "$O/UPS/Samples/""$Current_Patient_id/Copy_Number_Variation" ]; then
        mkdir -p "$O/UPS/Samples/$Current_Patient_id/Copy_Number_Variation"
    fi

    if [ ! -d "$O/UPS/Samples/""$Current_Patient_id/Copy_Number_Variation/CNVkit" ]; then
        mkdir -p "$O/UPS/Samples/$Current_Patient_id/Copy_Number_Variation/CNVkit"
    fi

    if [ ! -d "$O/UPS/Samples/""$Current_Patient_id/Plots" ]; then
        mkdir -p "$O/UPS/Samples/$Current_Patient_id/Plots"
    fi
done 

#################################### 
# Get BAMs for each Sample
####################################

for Current_Patient_id in "${unique_patient_ids[@]}"; do
    for line in "${SM_array[@]}"; do
        IFS=',' read -ra fields <<< "$line"
        # Check if the line contains the Patient_ID and "SRA" or "FASTq"
        
        if [[ "$line" == *"$Current_Patient_id"*"SRR"* ]]; then
            filePath=$(echo "$line" | awk '{print $11}')
            CurrentSample=$(echo "$line" | awk '{print $1}')
            cd "$O/UPS/Samples/$Current_Patient_id/FastQ"
            ReferenceGenome="$Support_Path/HG38.fasta"
            

            if [ ! -f "$CurrentSample""_pass_1.fastq.gz" ]; then
                echo "Extraction Needed for $CurrentSample"
                while [[ $(Count_Total_Jobs) -gt $Batch_Job_Limit ]]; do
                    sleep 30
                done
                sbatch "$UPS_Script_File""UPS_Fastq-dump.sh" "$filePath"
            fi
        fi

        if [[ "$line" == *"$Current_Patient_id"*"FastQ"* ]]; then
            filePath=$(echo "$line" | awk -F',' '{print $10}')
            cd "$O/UPS/Samples/$Current_Patient_id/FastQ"
            # Remove the _1.fastq.gz substring if it is still there
            if [[ $filePath == *_1.fastq.gz ]]; then
                filePath="${filePath%%_1.fastq.gz}"
            fi
            filePath="$(echo -e "${filePath}" | tr -d '[:space:]')"
            if [ ! -f "$filePath""_1.fastq.gz" ]; then
                cp $Input_Path$filePath""_1.fastq.gz .
                cp $Input_Path$filePath""_2.fastq.gz .
            fi
           
        fi
        

    

        if [[ "$line" == *"$Current_Patient_id"*"BAM"* ]]; then
            cd "$O/UPS/Samples/$Current_Patient_id/BAM"
            filePath=$(echo "$line" | awk '{print $10}')
        
            if [[ "${fields[4]}" == *"Yes"* ]]; then
                filename="$Current_Patient_id""_Tumour.bam"
                if [ ! -f "$filename" ]; then
                    echo "Moving $filename"
                    cp "$filePath" "$filename"
                fi
            fi
        
            if [[ "${fields[4]}" == *"No"* ]]; then
                filename="$Current_Patient_id""_Normal.bam"
                if [ ! -f "$filename" ]; then
                    echo "Moving $filename"
                    cp "$filePath" "$filename"
                fi
            fi
        fi
    done
done

Count_Job_Type "Fastq-dump"
while [ $Current_Job_Count -gt 0 ]; do
    sleep 30
done


for Current_Patient_id in "${unique_patient_ids[@]}"
do
    for line in "${SM_array[@]}"; do
        # Check if the line contains the Patient_ID and "SSR" or "FASTq"
        if [[ "$line" == *"$Current_Patient_id"*"SRR"* ]] || [[ $filePath == *"_1.fastq.gz"* ]]; then
            filePath=$(echo "$line" | awk '{print $10}')
            

            cd "$O/UPS/Samples/$Current_Patient_id/BAM"
            filePath="$(echo -e "${Current_Patient_id}" | tr -d '[:space:]')"
            File1="$Current_Patient_id""_1.fastq.gz"c
            File2="$Current_Patient_id""_2.fastq.gz"
            while [[ $(Count_Total_Jobs) -gt $Batch_Job_Limit ]]; do
                sleep 30
            done
            sbatch "$UPS_Script_File""UPS_Fastq_To_Bam.sh" "$Reference_Genome" "$File1" "$File2" "$Current_Sample"
        fi
    done
done 

Count_Job_Type "FastqToBam"
while [ $Current_Job_Count -gt 0 ]; do
    sleep 30
done

#################################### 
# Get SNV, SV, and CNVs for each Sample
####################################


for Current_Patient_id in "${unique_patient_ids[@]}"; do
    cd "$O/UPS/Samples/$Current_Patient_id/BAM"
    TumourFileFlag="Negtive"
    NormalFileFlag="Negtive"
    TumourFile=""
    NormalFile=""

    
    echo $Current_Patient_id
    for line in "${SM_array[@]}"; do
        IFS=',' read -ra fields <<< "$line"
        if [[ "$line" == *"$Current_Patient_id"* ]]; then
            if [[ "${fields[4]}" == *"Yes"* ]]; then
                TumourFileFlag="Present"
                TumourFile="$O/UPS/Samples/$Current_Patient_id/BAM/$Current_Patient_id""_Tumour.bam"
            elif [[ "${fields[4]}" == *"No"* ]]; then
                NormalFileFlag="Present"
                NormalFile="$O/UPS/Samples/$Current_Patient_id/BAM/$Current_Patient_id""_Normal.bam"
            fi
    
            if [[ "${fields[5]}" == *"WGS"* ]]; then
                WGSWESFlag="WGS"
            elif [[ "${fields[5]}" == *"WES"* ]]; then
                WGSWESFlag="WES"
            fi
        fi
    done

    #check if SNV have allready been called 
    if [ ! -f "$O/UPS/Samples/$Current_Patient_id/Single_Nucleotide_Variant/Mutect/Filtered_somatic.vcf.gz" ]; then
        while [[ $(Count_Total_Jobs) -gt $Batch_Job_Limit ]]; do
            sleep 30
        done
        cd "$O/UPS/Samples/$Current_Patient_id/Single_Nucleotide_Variant/Mutect/"
        sbatch "$UPS_Script_File""UPS_SNV_Mutect.sh" \
        -R $Reference_Genome -S $Current_Patient_id -ST $Sequencing_Type \
        -PTFlaf $TumourFileFlag -PTFile $TumourFile \
        -PNFlaf $NormalFileFlag -PNFile $NormalFile \
        -Support $Support_Path -O $O 
        echo "sbatch "$UPS_Script_File""UPS_SNV_Mutect.sh" \
        -R $Reference_Genome -S $Current_Patient_id -ST $Sequencing_Type \
        -PTFlaf $TumourFileFlag -PTFile $TumourFile \
        -PNFlaf $NormalFileFlag -PNFile $NormalFile \
        -Support $Support_Path -O $O "
    fi
    
    #check if SV have allready been called 
    if [ ! -f "$O/UPS/Samples/$Current_Patient_id/Structural_Variant/Manta/Manta.vcf" ]; then
        while [[ $(Count_Total_Jobs) -gt $Batch_Job_Limit ]]; do
            sleep 30
        done
        cd "$O/UPS/Samples/$Current_Patient_id/Structural_Variant/Manta/"
        sbatch "$UPS_Script_File""UPS_SV_Manta.sh" \
        -R "$Reference_Genome" -S "$Current_Patient_id" -ST "$Sequencing_Type" \
        -PTFlaf "$TumourFileFlag" -PTFile "$TumourFile" \
        -PNFlaf "$NormalFileFlag" -PNFile "$NormalFile" \
        -Support $Support_Path -O $O 
        echo "sbatch "$UPS_Script_File""UPS_SV_Manta.sh" \
        -R "$Reference_Genome" -S "$Current_Patient_id" -ST "$Sequencing_Type" \
        -PTFlaf "$TumourFileFlag" -PTFile "$TumourFile" \
        -PNFlaf "$NormalFileFlag" -PNFile "$NormalFile" \
        -Support $Support_Path -O $O squeue"
    fi

     #check if CNV have allready been called 
    if [ ! -f "$O/UPS/Samples/$Current_Patient_id/Copy_Number_Variation/CNVkit" ]; then
        while [[ $(Count_Total_Jobs) -gt $Batch_Job_Limit ]]; do
            sleep 30
        done
        cd "$O/UPS/Samples/$Current_Patient_id/Copy_Number_Variation/CNVkit"
        sbatch "$UPS_Script_File""UPS_CNV_CNVkit.sh" \
        -R "$Reference_Genome" -S "$Current_Patient_id" \
        -PTFlaf "$TumourFileFlag" -PTFile "$TumourFile" \
        -PNFlaf "$NormalFileFlag" -PNFile "$NormalFile" \
        -Support $Support_Path -O $O -WE $WGSWESFlag 
        
        echo "sbatch "$UPS_Script_File""UPS_CNV_CNVkit.sh" \
        -R "$Reference_Genome" -S "$Current_Patient_id" \
        -PTFlaf "$TumourFileFlag" -PTFile "$TumourFile" \
        -PNFlaf "$NormalFileFlag" -PNFile "$NormalFile" \
        -Support $Support_Path -O $O -WE $WGSWESFlag "
       
    fi    

    # Summerise result in gene centric tables 

    # Visualization 

    # Sample Report 

done 

# Summery tables

# Visualization

# Total Report


