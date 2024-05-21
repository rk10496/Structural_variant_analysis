#!/bin/bash

# This script assigns all reads in a single read group, marks duplicate reads, and recalibrates the base score.

# Author=Rahul_kumar/LeeOestereich lab


#SBATCH --job-name=Mapping_recalibration
#SBATCH -N 1
#SBATCH --cpus-per-task=64
#SBATCH -t 2-00:00
 
# Modules required
module purge
Picard=/ihome/crc/install/picard/2.18.12
GATK=/ihome/crc/install/gatk/gatk-4.5.0.0
module load gcc/8.2.0
 
#Create a new directory
designated_path="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample"
mkdir -p "$designated_path/5_mkdup"
 
#Input directory
input_dir1="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/4_mapped_file"
input_dir2="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/5_mkdup"

# Output directory
output_dir="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/5_mkdup"

# Reference genome fasta file
reference="/bgfs/alee/LO_LAB/Personal/Rahul/Software/Human_ref_GRCH38.p7/Human_genome_hg38_p7_chr_only.fna"

# Known sites for BQSR 
known_sites="/bgfs/alee/LO_LAB/Personal/Rahul/Software"

# AddOrReplaceReadGroups, assigns all the reads in a file to a single new read-group.

# Iterate over each sorted BAM file in the input directory
for sorted_bam in "$input_dir1"/*_sorted.bam; do

    # Extract sample name from the sorted BAM file name
    sample=$(basename "$sorted_bam" _sorted.bam)

    # Output RG BAM file for the sample
    rg_bam="${output_dir}/${sample}.rg.bam"
   
    # Run AddOrReplaceReadGroups tool
    java -jar $Picard/picard.jar AddOrReplaceReadGroups \
        I="$sorted_bam" \
        O="$rg_bam" \
        RGID=1 \
        RGLB=library \
        RGPL=illumina \
        RGPU=unit \
        RGSM="$sample"
    
    echo "Processed BAM file created: $rg_bam"
	
done

echo "Processing completed for all BAM files."

# MarkDuplicates by Picard, locates and tags duplicate reads in a BAM file.

# Iterate over each RG BAM file in the input directory
for rg_bam in "$input_dir2"/*_rg.bam; do

    # Extract sample name from the RG BAM file name
    sample=$(basename "$rg_bam" _rg.bam)

    # Output marked duplicates BAM file for the sample
	rg_mkdp_bam="${output_dir}/${sample}.rg.mkdp.bam"
  

    # Run Picard's MarkDuplicates tool
    java -jar $Picard/picard.jar MarkDuplicates \
        I="$rg_bam" \
        O="$rg_mkdp_bam" \
        M="$output_dir/${sample}rg_dup_metrics.txt" \
        REMOVE_DUPLICATES=false

    echo "Marked duplicates BAM file created: $rg_mkdp_bam"

done

echo "Marking duplicates completed for all samples."

# BuildBamIndex, used to generate a BAM index ".bai" file.

# Iterate over each rg_mkdp_bam file in the input directory
for rg_mkdp_bam in "$input_dir2"/*_rg.mkdp.bam; do

    # Extract sample name from the rg_mkdp_bam file name
    sample=$(basename "$rg_mkdp_bam" _rg.mkdp.bam)

    # Output rg mkdp bam bai file for the sample
     rg_mkdp_bam_bai="${output_dir}/${sample}.rg.mkdp.bam.bai"
   
    # Run BuildBamIndex utility
    java -jar $Picard/picard.jar BuildBamIndex \
        I="$rg_mkdp_bam" \
        O="$rg_mkdp_bam_bai" \
    
    echo "Index file created: $rg_mkdp_bam_bai"
		
done

echo "Indexing completed for all rg_mkdup_bam files."


# BaseRecalibrator used to generates recalibration table for Base Quality Score Recalibration (BQSR).
# ApplyBQSR to recalibrates the base qualities of the input reads based on the recalibration table.


# Iterate over each rg_mkdp_bam file in the input directory
for rg_mkdp_bam in "$input_dir"/*rg.mkdp.bam; do

    # Extract sample name from the rg_mkdp_bam file name
    sample=$(basename "$rg_mkdp_bam" rg.mkdp.bam)
    
    # Output recalibration table for the sample
    recal_table="${output_dir}/${sample}.recal.table"
    
    # Output recalibrated BAM file
    recal_bam="${output_dir}/${sample}.recal.bam"
    
    # Run BaseRecalibrator utility
    java -jar /ihome/crc/install/gatk/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar BaseRecalibrator \
        -I "$rg_mkdp_bam" \
        -R "$reference" \
        --known-sites "$known_sites"/Homo_sapiens_assembly38.known_indels.vcf.gz \
        --known-sites "$known_sites"/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
        -O "$recal_table"
    
    echo "Recalibration table created for $sample"
    
    # Apply BQSR utility
    java -jar /ihome/crc/install/gatk/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar ApplyBQSR \
        -I "$rg_mkdp_bam" \
        -R "$reference" \
        --bqsr-recal-file "$recal_table" \
        -O "$recal_bam"
    
    echo "Recalibrated BAM file created for $sample"
done

echo "Recalibrated BAM files generated for all rg_mkdup_bam files."

# BuildBamIndex, used to generate a BAM index ".bai" file.

# Iterate over each recal_bam file in the input directory
for recal_bam in "$input_dir2"/*.recal.bam; do

    # Extract sample name from the recal_bam file name
    sample=$(basename "$recal_bam" .recal.bam)

    # Output rg mkdp bam bai file for the sample
     recal_bam_bai="${output_dir}/${sample}.recal.bam.bai"
   
    # Run BuildBamIndex utility
    java -jar $Picard/picard.jar BuildBamIndex \
        I="$recal_bam" \
        O="$recal_bam_bai" \
    
    echo "Index file created: $recal_bam_bai"
		
done

echo "Indexing completed for all recal_bam files."
