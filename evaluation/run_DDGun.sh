#!/bin/bash

for i in $(seq $1 $2);
  do python3 ./ddgun_seq.py /synthetic_point_mutations/$i.fasta /synthetic_point_mutations/$i.muts > $i.ddgun_results;
done