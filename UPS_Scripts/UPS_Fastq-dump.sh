#!/bin/bash
#SBATCH --job-name=Fastq-dump
module load SRA
SSRPath="$1"
echo "$SSRPath"
fastq-dump --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip "$SSRPath"


