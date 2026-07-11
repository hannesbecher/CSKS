#!/bin/bash
# Sweep the PCR amplification-noise parameter at fixed depth/complexity and
# record the empirical duplication rate, analogous to run_duplication_analysis.sh
# but using 06delPipeCoalComplex.sh / makeReadsImperfectPCR.py.
ulimit -n 2048

CSKS_DIR="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
TMP_BASE="${CSKS_DIR}/tmp"
mkdir -p "$TMP_BASE"
OUTPUT="$CSKS_DIR/imperfect_pcr_results.csv"

# fixed simulation settings -- adjust depth/cplx here to match whatever
# combination you want to profile the amplification-noise effect at
depth=20
cplx=0.5

echo "pcr_dist,pcr_param,depth,cplx,r,dup_empirical" > "$OUTPUT"

run_one() {
  local dist="$1"
  local param="$2"

  echo "Running pcr_dist=$dist pcr_param=$param depth=$depth cplx=$cplx"

  bash "$CSKS_DIR/06delPipeCoalComplex.sh" \
    -g 1000000 \
    -t 0.01 \
    --depth "$depth" \
    --cplx "$cplx" \
    --pcr-dist "$dist" \
    --pcr-param "$param" \
    --pref "G1000000P2T0.01D${depth}C${cplx}PCR${dist}${param}" \
    -d "$TMP_BASE" \
    --keep 2>/dev/null

  reads=$(ls -t ${TMP_BASE}/csks.*/reads.fa 2>/dev/null | head -1)
  if [ -z "$reads" ]; then
    echo "WARNING: no reads found for pcr_dist=$dist pcr_param=$param"
    return
  fi

  total=$(grep -c "^>" "$reads")
  unique=$(grep -v "^>" "$reads" | sort | uniq | wc -l | tr -d ' ')

  r=$(python3 -c "print(round((1000000*${depth}*2/150)/(${cplx}*1990000), 4))")
  dup=$(python3 -c "print(round(1 - ${unique}/${total}, 4))")

  echo "$dist,$param,$depth,$cplx,$r,$dup" >> "$OUTPUT"
  echo "  r=$r | total=$total | unique=$unique | dup=$dup"

  rm -rf "$reads"
  rmdir "$(dirname "$reads")" 2>/dev/null
}

# baseline: no amplification noise (should match the "perfect PCR" duplication
# rate already recorded in duplication_rates.csv / duplication_results.csv for
# the same depth/cplx combination). run_one() isn't used here because
# --pcr-dist none must be called without --pcr-param, unlike geometric/lognormal.
echo "Running pcr_dist=none (baseline) depth=$depth cplx=$cplx"
bash "$CSKS_DIR/06delPipeCoalComplex.sh" \
  -g 1000000 -t 0.01 --depth "$depth" --cplx "$cplx" \
  --pcr-dist none \
  --pref "G1000000P2T0.01D${depth}C${cplx}PCRnone" \
  -d "$TMP_BASE" --keep 2>/dev/null
reads=$(ls -t ${TMP_BASE}/csks.*/reads.fa 2>/dev/null | head -1)
total=$(grep -c "^>" "$reads")
unique=$(grep -v "^>" "$reads" | sort | uniq | wc -l | tr -d ' ')
r=$(python3 -c "print(round((1000000*${depth}*2/150)/(${cplx}*1990000), 4))")
dup=$(python3 -c "print(round(1 - ${unique}/${total}, 4))")
echo "none,,$depth,$cplx,$r,$dup" >> "$OUTPUT"
echo "  r=$r | total=$total | unique=$unique | dup=$dup"
rm -rf "$reads"; rmdir "$(dirname "$reads")" 2>/dev/null

# geometric amplification, sweeping success probability p (p=1 -> no noise,
# decreasing p -> increasing amplification variance)
for p in 1.0 0.5 0.3 0.1 0.05 0.02; do
  run_one geometric "$p"
done

echo "Done. Results saved to $OUTPUT"
cat "$OUTPUT"
