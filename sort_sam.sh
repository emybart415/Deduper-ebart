#!/bin/bash

#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=interactive           #REQUIRED: which partition to use
#SBATCH --cpus-per-task=8                 #optional: number of cpus, default is 1
#SBATCH --mem=16GB                        #optional: amount of memory, default is 4GB


input_sam="/projects/bgmp/shared/deduper/C1_SE_uniqAlign.sam"
output_bam="/projects/bgmp/ebart/bioinfo/Bi624/Deduper-ebart/C1_SE_uniqAlign.bam"
sorted_bam="/projects/bgmp/ebart/bioinfo/Bi624/Deduper-ebart/C1_SE_uniqAlign_sorted.bam"
sorted_sam="/projects/bgmp/ebart/bioinfo/Bi624/Deduper-ebart/C1_SE_uniqAlign_sorted.sam"


samtools view -bS $input_sam -o $output_bam
samtools sort $output_bam -o $sorted_bam
samtools view -h $sorted_bam -o $sorted_sam

