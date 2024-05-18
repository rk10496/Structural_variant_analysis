#!/bin/bash

# This script is used to align the samples with the reference genome and check the quality of BAM files.

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
module load qualimap/2.2.2
module load  multiqc/1.19

# Create new directory
designated_path="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample"
mkdir -p "$designated_path/4_alignment_file"

designated_path1="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/4_alignment_file"
mkdir -p "$designated_path1/qualimap"

# Reference genome file
reference="/bgfs/alee/LO_LAB/Personal/Rahul/software/Human_ref_GRCH38.p7/Human_genome_hg38_p7_chr_only.fna"

# Directory containing FASTQ files for alignment
input_dir1="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/2_trimmed_file/paired"
input_dir2="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/4_alignment_file"
input_dir3="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/4_alignment_file/qualimap"

# Output directory for aligned SAM files
output_dir="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/4_alignment_file"
output_dir1="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/4_alignment_file/qualimap"
	
# List of sample names
samples=("TP19-M480_FOL6151A5_S12" "TP19-M483_FOL6151A4_S9" "TP19-M497_FOL6151A3_S7" "TP19-M774_FOL6151A1_S2" "TP19-M892_FOL6151A2_S4" )

# Iterate over each sample
for sample in "${samples[@]}"; do
	
    # input FASTQ files for the sample
    Forward_read="${input_dir1}/${sample}_R1_paired.fastq.gz"
    Reverse_read="${input_dir1}/${sample}_R2_paired.fastq.gz"

    # Output SAM file for the sample
    output_sam="${Output_dir}/${sample}.sam"

    # Perform alignment using bwa-mem
    echo "Aligning sample: $sample"
    bwa mem -t 64 $reference $forward_read $reverse_read > $output_sam

    echo "Alignment completed for $sample"
done

echo "All the samples were aligned successfully."

# samtools is a suite of programs for interacting with high-throughput sequencing data. 
# samtools view used to convert a sam file to a bam file
# samtools sort used to sort the sam/bam file by coordinates.
# samtools index used to create an index file for the sorted bam file.

# Iterate over each SAM file in the input directory
for sam_file in $input_dir2/*.sam; do
    # Extract sample name from the SAM file
    sample=$(basename "$sam_file" .sam)

    # Output BAM file
    bam_file="$output_dir/$sample.bam"

    # Convert SAM to BAM
    echo "Converting $sam_file to BAM format..."
    samtools view -bS "$sam_file" > "$bam_file"

    echo "BAM file created: $bam_file"
done

echo "BAM files generated for all the SAM files."

# Iterate over each BAM file in the input directory
for bam_file in $input_dir2/*.bam; do
    # Extract sample name from the BAM file
    sample=$(basename "$bam_file" .bam)

    # Output Sorted BAM file for the sample
    sorted_bam="${output_dir}/${sample}_sorted.bam"
	
	# Sort BAM file
    echo "Sorting $bam_file..."
    samtools sort "$bam_file" -o "$sorted_bam"

    echo "Sorted BAM file created: $sorted_bam"
done

echo "All samples were sorted successfully."

# Iterate over each Sorted BAM file in the input directory
for sorted_bam_file in $input_dir2/*_sorted.bam; do
    # Extract sample name from the sorted BAM file
    sample=$(basename "$sorted_bam_file" _sorted.bam)

    # Output bai file for the sample
    output_bai="${output_dir}/${sample}.bai"
	
	# Index BAM file
    echo "Indexing $sorted_bam_file..."
    samtools index "$sorted_bam_file"

    echo "Indexing bai file created: $Output_bai"
done

echo "All samples were converted sorted and indexed successfully."

#BamQC, used to perform quality control on BAM files

# Iterate over each BAM file in the input directory
for bam_file in "$input_dir2"/*sorted.bam; do
    # Extract sample name from the BAM file name
    sample=$(basename "$bam_file" _sorted.bam)

     # Create an output directory for the sample
    sample_output_dir="$output_dir1/$sample"
    mkdir -p "$sample_output_dir"

    # Run BamQC tool
    bamqc "$bam_file" --outdir "$sample_output_dir"

    echo "BamQC report generated for $sample"

    # Run Qualimap
    qualimap bamqc \
        -bam "$bam_file" \
        -outdir "$sample_output_dir" \
        -outfile "${sample}_qualimap_report.pdf" \
        -outformat pdf \
        --java-mem-size=4G

    echo "Qualimap report generated for $sample"
done

echo "QC completed for all samples"

# Load MultiQC module
module load multiqc/1.19

# Run MultiQC on the output directory
multiqc "$input_dir3" -o "$output_dir1"

echo "MultiQC analysis completed"
