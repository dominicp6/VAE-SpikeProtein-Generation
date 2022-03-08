#!/bin/bash

for i in $(seq $1 $2);
do python3 ./ddgun_seq.py ./synthetic_point_mutations/${i}_reference.fasta ./synthetic_point_mutations/$i.muts > ./synthetic_point_mutations/$i.txt;
done
