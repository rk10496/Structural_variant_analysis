#!/bin/bash

# This script is used to check the quality of high-throughput sequence data.

# Author=Rahul_kumar/LeeOestereich lab

#SBATCH --job-name=FastQC
#SBATCH -N 1
#SBATCH --cpus-per-task=64
#SBATCH -t 1-00:00

# Modules required
module purge
module load  fastqc/0.11.7
module load  multiqc/1.19

# Create new directory
designated_folder="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample"
mkdir -p "$designated_folder/1_fastqc_result"

# Input directory
input_dir="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample"

# Output directory 
output_dir="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/1_fastqc_result"

# Loop through all FASTQ files in the input directory
for fastq_file in $input_dir/*.fastq.gz; do
    # Extract the filename without extension
    sample=$(basename "$fastq_file" .fastq.gz)
    
    # Run FastQC
	fastqc -o "$output_dir" "$fastq_file" > "$output_dir/${sample}_fastqc.log"
    
    echo "FastQC analysis completed for $sample"
done

echo "All FastQC analyses completed"

# MultiQC is a reporting tool that aggregate results from bioinformatics analyses across many samples into a single report.

# Define input directory containing analysis results
input_dir1="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/1_fastqc_result"

# Define output directory where MultiQC report will be saved
output_dir1="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/1_fastqc_result"

# Run MultiQC on the input directory
multiqc "$input_dir1" -o "$output_dir1"

echo "MultiQC analysis completed"
