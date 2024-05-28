#!/bin/bash

# Author=Rahul_kumar/LeeOestereich lab

#SBATCH --job-name=Delly
#SBATCH -N 1
#SBATCH --cpus-per-task=64
#SBATCH -t 3-00:00


# Modules required
module purge
module load gcc/8.2.0
module load delly/1.0.3

# Create new directory
designated_path="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample"
mkdir -p "$designated_path/Delly"

# Reference genome file
reference="/bgfs/alee/LO_LAB/Personal/Rahul/Software/Human_ref_GRCH38.p7/Human_genome_hg38_p7_chr_only.fna"

# Input directory for bam files
input_dir="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/5_mkdup/"

# Output directory for aligned SAM files
output_dir="/bgfs/alee/LO_LAB/Personal/Rahul/Test_sample/Delly"

# Iterate over each BAM file in the input directory
for input_bam in ${input_dir}/*.recal.bam; do
    # Extract sample name from BAM file name
    sample=$(basename "${input_bam}" .recal.bam)

# Run Delly to detect structural variants
delly call -g ${reference} -o ${output_dir}/${sample}.vcf ${input_bam}
delly call -t DEL -g ${reference} -o ${output_dir}/${sample}_Del.vcf ${input_bam}
delly call -t DUP -g ${reference} -o ${output_dir}/${sample}_Dup.vcf ${input_bam}
delly call -t INV -g ${reference} -o ${output_dir}/${sample}_Inv.vcf ${input_bam}
delly call -t TRA -g ${reference} -o ${output_dir}/${sample}_Tra.vcf ${input_bam}
delly call -t INS -g ${reference} -o ${output_dir}/${sample}_Ins.vcf ${input_bam}

done

echo "Delly analysis completed for all samples"

