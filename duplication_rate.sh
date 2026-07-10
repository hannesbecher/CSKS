#!/bin/bash
ulimit -n 2048
bash 05delPipeCoalComplex.sh -g 1000000 -t 0.01 --depth 20 --cplx 0.5 \
  --pref G1000000P2T0.01D20C0.5 \
  -d ~/Desktop/The_University_of_Edinburgh/Sem2/Dissertation/Tetmer/CSKS_simulations/tmp \
  --keep

# Find reads file
reads=$(ls -t ~/Desktop/The_University_of_Edinburgh/Sem2/Dissertation/Tetmer/CSKS_simulations/tmp/csks.*/reads.fa | head -1)

echo "Using: $reads"

# Total reads
total=$(grep -c "^>" "$reads")

# Unique reads by sequence
unique=$(grep -v "^>" "$reads" | sort | uniq | wc -l)

# Duplication rate
echo "Total reads: $total"
echo "Unique reads: $unique"
python3 -c "print('Duplication rate:', round(1 - $unique/$total, 3))"
