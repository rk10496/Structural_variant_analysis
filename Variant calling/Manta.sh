#!/bin/bash

# This script is used for the Manta Structural Variant Caller

# Author=Rahul_kumar/LeeOestereich lab

#SBATCH --job-name=Manta
#SBATCH -N 1
#SBATCH --cpus-per-task=64
#SBATCH -t 1-00:00


# Modules required
module purge
module load gcc/8.2.0
module load manta/1.6.0

# Create new directory
designated_path="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample"
mkdir -p "$designated_path/Manta"

# Reference genome file
reference="/bgfs/alee/LO_LAB/Personal/Rahul/Software/Human_ref_GRCH38.p7/Human_genome_hg38_p7_chr_only.fna"

# Manta bin directory
manta_bin="/ihome/crc/install/manta/manta-1.6.0.centos6_x86_64/bin"

# Input directory for bam files
input_dir="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/5_mkdup/"

# Output directory 
output_dir="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/Manta"

# Iterate over each BAM file in the input directory
for input_bam in ${input_dir}/*.recal.bam; do
    # Extract sample name from BAM file name
    sample=$(basename "${input_bam}" .recal.bam)
	
	# Output directory for each sample
    sample_output_dir="$output_dir/$sample"
    mkdir -p "$sample_output_dir"
	
	 # Run the Manta configuration script
	 $manta_bin/configManta.py \
        --bam "$input_bam" \
        --referenceFasta "$reference" \
        --runDir "$sample_output_dir"
		
	 # Run Manta analysis
    $sample_output_dir/runWorkflow.py -m local -j 8

    echo "Manta analysis completed for sample: $sample"
done

echo "Manta analysis completed for all samples"	

