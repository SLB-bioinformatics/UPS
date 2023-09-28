#!/bin/bash
#SBATCH --job-name=Fastq-dump
module load SRA
filePath="$1"
echo "$filePath"
fastq-dump --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip "$$filePath"
