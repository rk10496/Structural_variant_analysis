#!/bin/bash

# This script is used to align the samples with the reference genome and check the qulaity of BAM files.

# Author=Rahul_kumar/LeeOestereich lab

# BWA is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome
# BWA-MEM, for high-quality queries and better performance, as it is faster and more accurate

#SBATCH --job-name=Alignment
#SBATCH -N 1
#SBATCH --cpus-per-task=64
#SBATCH -t 1-00:00

# Module required
module purge
module load gcc/8.2.0
module load bwa/0.7.17
module load samtools/1.9

# Create new directory
Designated_path="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample"
mkdir -p "$Designated_path/4_Mapped_file"

Designated_path1="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/4_Mapped_file"
mkdir -p "$Designated_path1/Qualimap"

# Reference genome file
Reference="/bgfs/alee/LO_LAB/Personal/Rahul/software/Human_ref_GRCH38.p7/GCF_000001405.33_GRCh38.p7_genomic.fna"

# Directory containing FASTQ files for mapping
Input_dir1="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/2_Trimmed_file/Paired"
Input_dir2="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/4_Mapped_file"
Input_dir3="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/4_Mapped_file/Qualimap"

# Output directory for aligned SAM files
Output_dir="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/4_Mapped_file"
Output_dir1="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/4_Mapped_file/Qualimap"
	
# List of sample names
Samples=("TP19-M480_FOL6151A5_S12" "TP19-M483_FOL6151A4_S9" "TP19-M497_FOL6151A3_S7" "TP19-M774_FOL6151A1_S2" "TP19-M892_FOL6151A2_S4" )


# Iterate over each sample
for Sample in "${Samples[@]}"; do
	
    # Input FASTQ files for the sample
    Forward_read="${Input_dir1}/${Sample}_R1_paired.fastq.gz"
    Reverse_read="${Input_dir1}/${Sample}_R2_paired.fastq.gz"

    # Output SAM file for the sample
    Output_sam="${Output_dir}/${Sample}.sam"

    # Perform mapping using bwa-mem
    echo "Aligning Sample: $Sample"
    bwa mem -t 64 $Reference $Forward_read $Reverse_read > $Output_sam

    echo "Alignment completed for $Sample"
done

echo "All the samples were mapped successfully."

# samtools is a suite of programs for interacting with high-throughput sequencing data. 
# samtools view used to convert a sam file to a bam file
# samtools sort used to sort the sam/bam file by coordinates.
# samtools index used to create an index file for the sorted bam file.

# Iterate over each SAM file in the input directory
for sam_file in $Input_dir2/*.sam; do
    # Extract sample name from the SAM file
    Sample=$(basename "$sam_file" .sam)

    # Output BAM file
    bam_file="$Output_dir/$Sample.bam"

    # Convert SAM to BAM
    echo "Converting $sam_file to BAM format..."
    samtools view -bS "$sam_file" > "$bam_file"

    echo "BAM file created: $bam_file"
done

echo "BAM files generated for all the SAM files."

# Iterate over each BAM file in the input directory
for bam_file in $Input_dir2/*.bam; do
    # Extract sample name from the BAM file
    Sample=$(basename "$bam_file" .bam)

    # Output Sorted BAM file for the sample
    sorted_bam="${Output_dir}/${Sample}_sorted.bam"
	
	# Sort BAM file
    echo "Sorting $bam_file..."
    samtools sort "$bam_file" -o "$sorted_bam"

    echo "Sorted BAM file created: $sorted_bam"
done

echo "All samples were sorted successfully."

# Iterate over each Sorted BAM file in the input directory
for sorted_bam_file in $Input_dir2/*_sorted.bam; do
    # Extract sample name from the sorted BAM file
    Sample=$(basename "$sorted_bam_file" _sorted.bam)

    # Output bai file for the sample
    Output_bai="${Output_dir}/${Sample}.bai"
	
	# Index BAM file
    echo "Indexing $sorted_bam_file..."
    samtools index "$sorted_bam_file"

    echo "Indexing bai file created: $Output_bai"
done

echo "All samples were converted sorted and indexed successfully."

#BamQC, used to perform quality control on BAM files

# Iterate over each BAM file in the input directory
for bam_file in "$Input_dir2"/*sorted.bam; do
    # Extract sample name from the BAM file name
    Sample=$(basename "$bam_file" _sorted.bam)

    # Create an output directory for the sample
    Sample_output_dir="$Output_dir1/$Sample"
    mkdir -p "$Sample_output_dir"

    # Run BamQC tool
    bamqc "$bam_file" --outdir "$Sample_output_dir"

    echo "BamQC report generated for $Sample"

    # Run Qualimap
    qualimap bamqc \
        -bam "$bam_file" \
        -outdir "$Sample_output_dir" \
        -outfile "${Sample}_qualimap_report.pdf" \
        -outformat pdf \
        --java-mem-size=4G

    echo "Qualimap report generated for $Sample"
done

echo "QC completed for all samples"

# Load MultiQC module
module load multiqc/1.19

# Run MultiQC on the output directory
multiqc "$Input_dir3" -o "$Output_dir1"

echo "MultiQC analysis completed"
