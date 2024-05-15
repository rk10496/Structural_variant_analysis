#!/bin/bash

# BWA is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome
# BWA-MEM, for high-quality queries and better performance, as it is faster and more accurate

#Author=Rahul_kumar

#SBATCH --job-name=Mapping
#SBATCH -N 1
#SBATCH --cpus-per-task=64
#SBATCH -t 1-00:00

# Module required
module purge
module load gcc/8.2.0
module load bwa/0.7.17
module load samtools/1.9

#Create new directory
Designated_path="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample"
mkdir -p "$Designated_path/4_Mapped_file"

# Reference genome file
Reference="/bgfs/alee/LO_LAB/Personal/Rahul/software/Human_ref_GRCH38.p7/GCF_000001405.33_GRCh38.p7_genomic.fna"

# Directory containing FASTQ files for mapping
Input_dir1="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/2_Trimmed_file/Paired"
Input_dir2="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/4_Mapped_file"

# Output directory for aligned SAM files
Output_dir="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/4_Mapped_file"
	
# List of sample names
Samples=("TP19-M480_FOL6151A5_S12" "TP19-M483_FOL6151A4_S9" "TP19-M497_FOL6151A3_S7" "TP19-M774_FOL6151A1_S2" "TP19-M892_FOL6151A2_S4" )

# Iterate over each sample
for sample in "${Samples[@]}"; do
	
    # Input FASTQ files for the sample
    Forward_read="${Input_dir1}/${sample}_R1_paired.fastq.gz"
    Reverse_read="${Input_dir1}/${sample}_R2_paired.fastq.gz"

    # Output SAM file for the sample
    Output_sam="${Output_dir}/${sample}.sam"

    # Perform mapping using bwa-mem
    echo "Aligning sample: $sample"
    bwa mem -t 64 $Reference $Forward_read $Reverse_read > $Output_sam

    echo "Mapping completed for sample: $sample"
done

echo "All the samples were mapped successfully."

# samtools is a suite of programs for interacting with high-throughput sequencing data. 
# samtools view used to convert a sam file to a bam file
# samtools sort used to sort the sam/bam file by coordinates.
# samtools index used to create an index file for the sorted bam file.

# Iterate over each SAM file in the input directory
for sam_file in $Input_dir2/*.sam; do
    # Extract sample name from the SAM file
    sample=$(basename "$sam_file" .sam)

    # Output BAM file
    bam_file="$Output_dir/$sample.bam"

    # Convert SAM to BAM
    echo "Converting $sam_file to BAM format..."
    samtools view -bS "$sam_file" > "$bam_file"

    echo "BAM file created: $bam_file"
done

echo "BAM files generated for all the SAM files."

# Iterate over each BAM file in the input directory
for bam_file in $Input_dir2/*.bam; do
    # Extract sample name from the BAM file
    sample=$(basename "$bam_file" .bam)

    # Output Sorted BAM file for the sample
    sorted_bam="${Output_dir}/${sample}_sorted.bam"
	
	# Sort BAM file
    echo "Sorting $bam_file..."
    samtools sort "$bam_file" -o "$sorted_bam"

    echo "Sorted BAM file created: $sorted_bam"
done

echo "All samples were sorted successfully."

# Iterate over each Sorted BAM file in the input directory
for sorted_bam_file in $Input_dir2/*_sorted.bam; do
    # Extract sample name from the sorted BAM file
    sample=$(basename "$sorted_bam_file" _sorted.bam)

    # Output bai file for the sample
    Output_bai="${Output_dir}/${sample}.bai"
	
	# Index BAM file
    echo "Indexing $sorted_bam_file..."
    samtools index "$sorted_bam_file"

    echo "Indexing bai file created: $Output_bai"
done

echo "All samples were converted sorted and indexed successfully.."
