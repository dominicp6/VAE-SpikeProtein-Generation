Reduced-size databases can be constructed from an original, large fasta file using the awk command as follows:

awk "(NR-1) % 100 == 0 || (NR-1) % 100 == 1" spikeprot0112.fasta > 1_in_50_spikeprot0112.fasta
awk "(NR-1) % 1000 == 0 || (NR-1) % 1000 == 1" spikeprot0112.fasta > 1_in_500_spikeprot0112.fasta
awk "(NR-1) % 10000 == 0 || (NR-1) % 10000 == 1" spikeprot0112.fasta > 1_in_5000_spikeprot0112.fasta