#!/bin/bash
ulimit -n 2048

CSKS_DIR="/Users/rahuldey/Desktop/The_University_of_Edinburgh/Sem2/Dissertation/Tetmer/CSKS_simulations/CSKS"
TMP_BASE="/Users/rahuldey/Desktop/The_University_of_Edinburgh/Sem2/Dissertation/Tetmer/CSKS_simulations/tmp"
OUTPUT="$CSKS_DIR/duplication_rates.csv"

echo "depth,cplx,r,dup_empirical" > "$OUTPUT"

for depth in 5 10 20 30; do
  for cplx in 0.1 0.2 0.4 0.6 0.8 1.0; do
    echo "Running depth=$depth cplx=$cplx"

    # Run simulation keeping temp files
    bash "$CSKS_DIR/05delPipeCoalComplex.sh" \
      -g 1000000 \
      -t 0.01 \
      --depth $depth \
      --cplx $cplx \
      --pref "G1000000P2T0.01D${depth}C${cplx}" \
      -d "$TMP_BASE" \
      --keep 2>/dev/null

    # Find the most recent temp directory
    reads=$(ls -t ${TMP_BASE}/csks.*/reads.fa 2>/dev/null | head -1)

    if [ -z "$reads" ]; then
      echo "WARNING: no reads found for depth=$depth cplx=$cplx"
      continue
    fi

    # Count total and unique reads
    total=$(grep -c "^>" "$reads")
    unique=$(grep -v "^>" "$reads" | sort | uniq | wc -l | tr -d ' ')

    # Calculate duplication rate and r
    r=$(python3 -c "print(round((1000000*${depth}*2/150)/(${cplx}*1990000), 4))")
    dup=$(python3 -c "print(round(1 - ${unique}/${total}, 4))")

    echo "$depth,$cplx,$r,$dup" >> "$OUTPUT"
    echo "  r=$r | total=$total | unique=$unique | dup=$dup"

    # Clean up this temp directory to save space
    rm -rf "$reads"
    rmdir "$(dirname $reads)" 2>/dev/null

  done
done

echo "Done. Results saved to $OUTPUT"
cat "$OUTPUT"
