#!/bin/bash
ulimit -n 2048

CSKS_DIR="/Users/rahuldey/Desktop/The_University_of_Edinburgh/Sem2/Dissertation/Tetmer/CSKS_simulations/CSKS"
TMP_BASE="/Users/rahuldey/Desktop/The_University_of_Edinburgh/Sem2/Dissertation/Tetmer/CSKS_simulations/tmp"
OUTPUT="$CSKS_DIR/duplication_results.csv"

DEPTHS="5 10 20 30"
COMPLEXITIES="0.01 0.02 0.04 0.08 0.16 0.32"
NREP=5
GENOME=1000000
THETA=0.01
PLOIDY=1
PDEL=0.0
RLEN=150

echo "depth,cplx,replicate,r,dup_empirical" > "$OUTPUT"

for depth in $DEPTHS; do
  for cplx in $COMPLEXITIES; do
    for rep in $(seq 1 $NREP); do

      pref="G${GENOME}P${PLOIDY}T${THETA}D${depth}C${cplx}R${rep}"
      echo "Running depth=$depth cplx=$cplx rep=$rep"

      # Run simulation keeping temp files
      bash "$CSKS_DIR/05delPipeCoalComplex.sh" \
        -g $GENOME \
        -t $THETA \
        -p $PLOIDY \
        --pdel $PDEL \
        --depth $depth \
        --cplx $cplx \
        --pref "$pref" \
        -d "$TMP_BASE" \
        --keep > /dev/null 2>&1

      # Get temp directory from log file
      logfile="$CSKS_DIR/${pref}.histXX.withdel${PDEL}.log"
      tmpdir=$(grep "^Temporary directory:" "$logfile" | awk '{print $3}')
      reads="$tmpdir/reads.fa"

      if [ ! -f "$reads" ]; then
        echo "WARNING: reads not found at $reads"
        echo "$depth,$cplx,$rep,NA,NA" >> "$OUTPUT"
        continue
      fi

      # Count total and unique reads
      total=$(grep -c "^>" "$reads")
      unique=$(grep -v "^>" "$reads" | sort | uniq | wc -l | tr -d ' ')

      # Calculate r and empirical duplication rate
      r=$(python3 -c "
n_reads  = $GENOME * $depth * $PLOIDY / $RLEN
n_unique = $cplx * $GENOME * $PLOIDY
print(round(n_reads / n_unique, 4))
")
      dup=$(python3 -c "print(round(1 - $unique/$total, 4))")

      echo "$depth,$cplx,$rep,$r,$dup" >> "$OUTPUT"
      echo "  r=$r | total=$total | unique=$unique | dup=$dup"

      # Clean up temp directory
      rm -rf "$tmpdir"

    done
  done
done

echo ""
echo "Done. Results saved to $OUTPUT"
cat "$OUTPUT"
