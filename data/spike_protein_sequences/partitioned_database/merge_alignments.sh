#!/bin/bash

path_to_muscle=$1
file_prefix_name=$2
start_number=$3
end_number=$4

counter=$start_number
echo Aligning sequences
while [ $counter -le $end_number ]
do
  echo Merging file number $counter
  $path_to_muscle -profile -in1 ${file_prefix_name}_${start_number},${end_number}.fasta -in2 ${file_prefix_name}_${counter}.fasta -out ${file_prefix_name}_${start_number},${end_number}.fasta
  ((counter++))
done
echo All done

#/home/dominic/miniconda3/pkgs/muscle-3.8.1551-h7d875b9_6/bin/muscle -profile -in1 spikeprot_0,25.afa -in2 spikeprot_26,45.afa -out spikeprot_0,45.afa
#/home/dominic/miniconda3/pkgs/muscle-3.8.1551-h7d875b9_6/bin/muscle -profile -in1 spikeprot_66,105.afa -in2 spikeprot_106,111.afa -out spikeprot_66,111.afa
#/home/dominic/miniconda3/pkgs/muscle-3.8.1551-h7d875b9_6/bin/muscle -profile -in1 spikeprot_0,65.afa -in2 spikeprot_66,111.afa -out spikeprot_all_aligned.afa