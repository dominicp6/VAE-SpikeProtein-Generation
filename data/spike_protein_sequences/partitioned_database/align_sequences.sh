#!/bin/bash

path_to_muscle=$1
file_prefix_name=$2
start_number=$3
end_number=$4

counter=$start_number
echo Aligning sequences
while [ $counter -le $end_number ]
do
  echo Starting file number $counter
  $path_to_muscle -in ${file_prefix_name}_${counter}.fasta -out ${file_prefix_name}_${counter}.afa
  ((counter++))
done
echo All done