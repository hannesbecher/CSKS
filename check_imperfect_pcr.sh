#!/bin/bash
# Self-contained smoke test for the imperfect-PCR extension.
# Skips KMC/histograms entirely -- just runs the genome + read simulation
# steps on a tiny genome and checks that the duplication rate responds to
# the PCR amplification-noise parameter the way it should.
#
# USAGE: run from the CSKS project root (the directory containing code/):
#   bash check_imperfect_pcr.sh
# Paste the full output back to Claude for a pass/fail check.

set -e
script_dir=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
CODE_DIR="$script_dir/code"

gs=5000
ploy=2
theta=0.01
pDel=0.0
cplx=0.5
depth=15

TMP_BASE=$(mktemp -d /tmp/pcrcheck.XXXXXX)
echo "Working directory: $TMP_BASE"
echo ""

RESULTS_FILE="$TMP_BASE/results.csv"
echo "dist,param,mean_amp,cv_amp,total_reads,unique_reads,dup_rate" > "$RESULTS_FILE"

run_case() {
  local dist="$1"
  local param="$2"
  local td="$TMP_BASE/${dist}_${param}"
  mkdir -p "$td"

  python3 "$CODE_DIR/makeGenomesDel.py" "$gs" "$ploy" "$theta" "$td" "$pDel" > "$td/genome.log" 2>&1

  local pcr_args=(--pcr-dist "$dist")
  if [ "$dist" != "none" ]; then
    pcr_args+=(--pcr-param "$param")
  fi

  python3 "$CODE_DIR/makeReadsImperfectPCR.py" "$td" "$cplx" --depth "$depth" "${pcr_args[@]}" > "$td/reads.log" 2>&1

  local mean_amp cv_amp
  mean_amp=$(grep "Amplification factor summary" "$td/reads.log" | sed -E 's/.*mean=([0-9.eE+-]+).*/\1/')
  cv_amp=$(grep "Amplification factor CV" "$td/reads.log" | sed -E 's/.*CV: ([0-9.eE+-]+).*/\1/')
  [ -z "$mean_amp" ] && mean_amp="NA"
  [ -z "$cv_amp" ] && cv_amp="NA"

  local total unique
  total=$(grep -c "^>" "$td/reads.fa")
  unique=$(grep -v "^>" "$td/reads.fa" | sort | uniq | wc -l | tr -d ' ')
  local dup
  dup=$(python3 -c "print(round(1 - $unique/$total, 4))")

  echo "$dist,$param,$mean_amp,$cv_amp,$total,$unique,$dup" >> "$RESULTS_FILE"
  echo "  $dist (param=$param): mean_amp=$mean_amp cv=$cv_amp total=$total unique=$unique dup_rate=$dup"

  rm -rf "$td"
}

echo "Running baseline (no PCR noise)..."
run_case none NA

echo ""
echo "Running geometric sweep (p: 1.0 -> 0.05)..."
for p in 1.0 0.5 0.2 0.05; do
  run_case geometric "$p"
done

echo ""
echo "Running lognormal sweep (sigma: 0.0 -> 2.0)..."
for s in 0.0 0.5 1.0 2.0; do
  run_case lognormal "$s"
done

echo ""
echo "########################################"
echo "RESULTS TABLE"
echo "########################################"
if command -v column >/dev/null 2>&1; then
  column -s, -t "$RESULTS_FILE"
else
  awk -F, '{printf "%-10s %-8s %-10s %-8s %-12s %-13s %-9s\n", $1,$2,$3,$4,$5,$6,$7}' "$RESULTS_FILE"
fi

echo ""
echo "########################################"
echo "AUTOMATIC CHECKS"
echo "########################################"
python3 - "$RESULTS_FILE" <<'PYEOF'
import csv, sys

path = sys.argv[1]
rows = list(csv.DictReader(open(path)))

def find(dist, param):
    for r in rows:
        if r["dist"] == dist and (dist == "none" or float(r["param"]) == float(param)):
            return r
    return None

checks = []

base = find("none", None)
geom1 = find("geometric", 1.0)
logn0 = find("lognormal", 0.0)

def close(a, b, tol=0.08):
    return abs(float(a) - float(b)) <= tol

checks.append((
    "none baseline dup rate close to geometric p=1.0 (no noise)",
    close(base["dup_rate"], geom1["dup_rate"]),
    f"none={base['dup_rate']} geom(p=1.0)={geom1['dup_rate']}"
))
checks.append((
    "none baseline dup rate close to lognormal sigma=0.0 (no noise)",
    close(base["dup_rate"], logn0["dup_rate"]),
    f"none={base['dup_rate']} lognorm(sigma=0)={logn0['dup_rate']}"
))

geom_rows = sorted([r for r in rows if r["dist"] == "geometric"], key=lambda r: -float(r["param"]))
geom_dups = [float(r["dup_rate"]) for r in geom_rows]
geom_monotonic = all(geom_dups[i] <= geom_dups[i+1] + 1e-9 for i in range(len(geom_dups)-1))
checks.append((
    "geometric: dup rate increases as p decreases (more noise)",
    geom_monotonic,
    "p/dup: " + ", ".join(f"{r['param']}={r['dup_rate']}" for r in geom_rows)
))

logn_rows = sorted([r for r in rows if r["dist"] == "lognormal"], key=lambda r: float(r["param"]))
logn_dups = [float(r["dup_rate"]) for r in logn_rows]
logn_monotonic = all(logn_dups[i] <= logn_dups[i+1] + 1e-9 for i in range(len(logn_dups)-1))
checks.append((
    "lognormal: dup rate increases as sigma increases (more noise)",
    logn_monotonic,
    "sigma/dup: " + ", ".join(f"{r['param']}={r['dup_rate']}" for r in logn_rows)
))

all_pass = True
for name, ok, detail in checks:
    status = "PASS" if ok else "FAIL"
    if not ok:
        all_pass = False
    print(f"[{status}] {name}")
    print(f"       {detail}")

print("")
print("OVERALL:", "ALL CHECKS PASSED" if all_pass else "SOME CHECKS FAILED")
PYEOF

rm -rf "$TMP_BASE"
