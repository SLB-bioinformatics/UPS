
Reference_Genome="${1}"
File1="${2}"
File2="${3}"
Current_Sample="${4}"

echo "$Reference_Genome"
echo "$File1"
echo "$File2"
echo "$Current_Sample"

# Check if all required arguments are provided
if [ $# -ne 4 ]; then
  echo "Usage: $0 <Reference_Genome> <File1> <File2> <Current_Sample>"
  exit 1        
fi

fastqc -o . $File1
fastqc -o . $File2

trimmomatic-0.40.jar PE -threads 8 \
    $File1 $File2 \
    $(echo "$File1" | sed 's/\.fastq.gz$/_paired.fastq.gz/') $(echo "$File1" | sed 's/\.fastq.gz$/_unpaired.fastq.gz/') \
    $(echo "$File2" | sed 's/\.fastq.gz$/_paired.fastq.gz/') $(echo "$File2" | sed 's/\.fastq.gz$/_unpaired.fastq.gz/') \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

File1=$(echo "$File1" | sed 's/\.fastq.gz$/_paired.fastq.gz/')
File2=$(echo "$File2" | sed 's/\.fastq.gz$/_paired.fastq.gz/')

# Aligning the input files
if [ ! -f "${Current_Sample}_1.sai" ]; then
  bwa aln "${Reference_Genome}" "${File1}" > "${Current_Sample}_1.sai"
fi

if [ ! -f "${Current_Sample}_2.sai" ]; then
  bwa aln "${Reference_Genome}" "${File2}" > "${Current_Sample}_2.sai"
fi

# Generating a SAM file by pairing the aligned sequences from both input files with the reference genome
bwa sampe "${Reference_Genome}" "${Current_Sample}_1.sai" "${Current_Sample}_2.sai" "${File1}" "${File2}" > "${Current_Sample}.sam"

# Convert SAM to BAM using Samtools
samtools view -S -b "${Current_Sample}.sam" > "${Current_Sample}.bam"

#sort the BAM file 
samtools sort "${Current_Sample}.bam" -o "${Current_Sample}_sorted.bam"
mv "${Current_Sample}_sorted.bam" "${Current_Sample}.bam"