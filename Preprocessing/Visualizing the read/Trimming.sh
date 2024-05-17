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
module load  fastqc/0.11.7
module load  multiqc/1.19

# Create new directory
Designated_path="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample"
mkdir -p "$Designated_path/2_Trimmed_file"
mkdir -p "$Designated_path/3_Trimmed_Fastqc_result"

# Array of multiple designated folders
Designated_folder1="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/2_Trimmed_file/"
Designated_folder2="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/3_Trimmed_Fastqc_result/"

# Array of multiple directories
Directory_name=("Paired" "Unpaired")

# Iterate over each designated folder
for folder in "${Designated_folder[@]}"; do
    # Iterate over each directory name
    for dir_name in "${Directory_name[@]}"; do
        # Create the directory in the current designated folder
        mkdir -p "$folder/$dir_name"
    done

# Input directory
Input_dir1="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/Raw_file/"
Input_dir2="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/2_Trimmed_file/Paired"
Input_dir3="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/3_Trimmed_Fastqc_result/Paired"

# Output directory 
Output_dir_paired="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/2_Trimmed_file/Paired"
Output_dir_Unpaired="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/2_Trimmed_file/Unpaired/"
Output_dir="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/3_Trimmed_Fastqc_result/Paired"

# Define adapter file
adapter_file="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/Trimmed_file/adapter.fa"

# Loop through all paired-end FASTQ files in the input directory
for forward_read in $Input_dir1/*_R1_001.fastq.gz; do
    # Extract the filename without extension
    filename=$(basename "$forward_read" _R1_001.fastq.gz)
	
    # Extract the filename without extension
    reverse_read="${filename}_R2_001.fastq.gz"   
	
    # Run Trimmomatic for paired-end reads
    module load trimmomatic/0.38 
	java -jar /ihome/crc/install/trimmomatic/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 64 -phred33 \
        $forward_read $reverse_read \
        "${Output_dir_paired}/${filename}_R1_paired.fastq.gz" "${Output_dir_Unpaired}/${filename}_R1_unpaired.fastq.gz" \
        "${Output_dir_paired}/${filename}_R2_paired.fastq.gz" "${Output_dir_Unpaired}/${filename}_R2_unpaired.fastq.gz" \
        ILLUMINACLIP:$adapter_file:2:30:10 LEADING:0 TRAILING:0 SLIDINGWINDOW:4:15 MINLEN:75
    
    echo "Trimming completed for $filename"
done

# FastQC,a tool used to check the quality of high-throughput sequence data
# MultiQC is a reporting tool that aggregate results from bioinformatics analyses across many samples into a single report.


# Loop through all FASTQ files in the input directory
for fastq_file in $Input_dir2/*paired.fastq.gz; do
    # Extract the filename without extension
    filename=$(basename "$fastq_file" paired.fastq.gz)
    
    # Run FastQC
    module load  fastqc/0.11.7
	fastqc -t 64 -f fastq -o $Output_dir $fastq_file
	
    echo "FastQC analysis completed for $filename"
done

echo "Analyses completed for all samples"

# Run MultiQC on the input directory
multiqc "$Input_dir3" -o "$Output_dir"

echo "MultiQC analysis completed"

done

echo "Trimming and Trimmed_FastQC completed for all samples"
