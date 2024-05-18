#!/bin/bash

# This script is used to trim low quality reads and perform QC analysis on trimmed read

# Author:Rahul_Kumar/LeeOestereich lab

# Before running trimmomatic, have a look at MultiQC Report
# Trimmomatic,trimming tool to remove adapter and low quality reads for Illumina NGS data.

#SBATCH --job-name=Trimming
#SBATCH -N 1
#SBATCH --cpus-per-task=64
#SBATCH -t 1-00:00

# Modules required
module purge
module load trimmomatic/0.38
module load  fastqc/0.11.7
module load  multiqc/1.19

# Create new directory
designated_path="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample"
mkdir -p "$designated_path/2_trimmed_file"
mkdir -p "$designated_path/3_trimmed_fastqc_result"

# Array of multiple designated folders
designated_folder1="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/2_trimmed_file/"
designated_folder2="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/3_trimmed_fastqc_result/"

# Array of multiple directories
directory_name=("paired" "unpaired")

# Iterate over each designated folder
for folder in "${designated_folder[@]}"; do
    # Iterate over each directory name
    for dir_name in "${directory_name[@]}"; do
        # Create the directory in the current designated folder
        mkdir -p "$folder/$dir_name"
    done

# Input directory
input_dir1="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/Raw_file/"
input_dir2="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/2_trimmed_file/paired"
input_dir3="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/3_trimmed_fastqc_result/paired"

# Output directory 
output_dir_paired="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/2_trimmed_file/paired"
output_dir_unpaired="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/2_trimmed_file/unpaired/"
output_dir="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/3_trimmed_fastqc_result/paired"

# Define adapter file
adapter_file="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/2_trimmed_file/adapter.fa"

# Loop through all paired-end FASTQ files in the input directory
for forward_read in $input_dir1/*_R1_001.fastq.gz; do
    # Extract the sample without extension
    sample=$(basename "$forward_read" _R1_001.fastq.gz)
	
    # Extract the sample without extension
    reverse_read="${sample}_R2_001.fastq.gz"   
	
    # Run Trimmomatic for paired-end reads 
	java -jar /ihome/crc/install/trimmomatic/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 64 -phred33 \
        $forward_read $reverse_read \
        "${output_dir_paired}/${sample}_R1_paired.fastq.gz" "${output_dir_unpaired}/${sample}_R1_unpaired.fastq.gz" \
        "${output_dir_paired}/${sample}_R2_paired.fastq.gz" "${output_dir_unpaired}/${sample}_R2_unpaired.fastq.gz" \
        ILLUMINACLIP:$adapter_file:2:30:10 LEADING:0 TRAILING:0 SLIDINGWINDOW:4:15 MINLEN:75
    
    echo "Trimming completed for $sample"
done

# FastQC,a tool used to check the quality of high-throughput sequence data
# MultiQC is a reporting tool that aggregate results from bioinformatics analyses across many samples into a single report.


# Loop through all FASTQ files in the input directory
for fastq_file in $input_dir2/*paired.fastq.gz; do
    # Extract the sample without extension
    sample=$(basename "$fastq_file" paired.fastq.gz)
    
    # Run FastQC
    module load  fastqc/0.11.7
	fastqc -t 64 -f fastq -o $output_dir $fastq_file
	
    echo "FastQC analysis completed for $sample"
done

echo "Analyses completed for all samples"

# Run MultiQC on the input directory
multiqc "$input_dir3" -o "$output_dir"

echo "MultiQC analysis completed"

done

echo "Trimming and Trimmed_FastQC completed for all samples"
