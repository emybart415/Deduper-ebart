#!/bin/bash

#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=interactive           #REQUIRED: which partition to use
#SBATCH --cpus-per-task=8                 #optional: number of cpus, default is 1
#SBATCH --mem=16GB                        #optional: amount of memory, default is 4GB
#SBATCH --error=deduper_%j.err

#input files
input_sam="/projects/bgmp/ebart/bioinfo/Bi624/Deduper-ebart/C1_SE_uniqAlign_sorted.sam"
umi="/projects/bgmp/ebart/bioinfo/Bi624/Deduper-ebart/STL96.txt"

#output directories
deduped_sam="/projects/bgmp/ebart/bioinfo/Bi624/Deduper-ebart/C1_SE_uniqAlign_sorted_deduped.sam"


/usr/bin/time -v python ./bartlett_deduper.py -u $umi -f $input_sam -o $deduped_sam
