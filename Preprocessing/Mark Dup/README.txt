We used Picard's (https://github.com/broadinstitute/picard) AddOrReplaceReadGroups, MarkDuplicates, and BuildBamIndex  to assign all the reads in a file to a single new read group, to locate and tag duplicate reads in a BAM file and indexing mkdup bam file respectively.
Followed by base recalibration with GATK BaseRecalibrator and Apply BQSR (https://github.com/broadinstitute/gatk).
