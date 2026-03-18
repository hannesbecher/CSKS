# activate conda environment (run 00setUpEnv.sh to set up)
#conda activate CSKS

usage() {
    cat <<EOF
Usage: $0 [-g genome_size] [-p ploidy] [-t theta] [--pdel proportion_deleted] [-d temp_base] [--keep] [-h]

Run the coalescent deletion pipeline and write k-mer histograms for multiple k values.

Options:
  -g, --gs       Chromosome size in nt (default: 10000000)
  -p, --ploy     Ploidy level (default: 2)
  -t, --theta    Population-scaled mutation rate (default: 0.1)
  --pdel         Proportion deleted (default: 0.01)
  -d, --tmpdir   Base path for the temporary directory (default: /tmp)
  --keep         Keep the temporary directory for debugging
  -h, --help     Show this help message and exit
EOF
}

# simulation parameters
gs=10000000 # Size in nt of each chromosome. There are 10 chromosomes, so the overall haploid genome size is 'gs'*10!
ploy=2 # ploidy level
theta=0.1 # population-scaled mutation rate (= heterozygosity in a random-mating population)
pDel=0.01 # proportion deleted
tmp_base=/tmp # base path for temporary directory
keep_tmp=0 # keep temporary directory for debugging

while [ "$#" -gt 0 ]; do
    case "$1" in
        -g|--gs)
            if [ -z "$2" ]; then
                echo "Missing value for $1" >&2
                usage >&2
                exit 1
            fi
            gs="$2"
            shift 2
            ;;
        -p|--ploy)
            if [ -z "$2" ]; then
                echo "Missing value for $1" >&2
                usage >&2
                exit 1
            fi
            ploy="$2"
            shift 2
            ;;
        -t|--theta)
            if [ -z "$2" ]; then
                echo "Missing value for $1" >&2
                usage >&2
                exit 1
            fi
            theta="$2"
            shift 2
            ;;
        --pdel)
            if [ -z "$2" ]; then
                echo "Missing value for $1" >&2
                usage >&2
                exit 1
            fi
            pDel="$2"
            shift 2
            ;;
        -d|--tmpdir)
            if [ -z "$2" ]; then
                echo "Missing value for $1" >&2
                usage >&2
                exit 1
            fi
            tmp_base="$2"
            shift 2
            ;;
        --keep)
            keep_tmp=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

output_prefix="G${gs}P${ploy}T${theta}"
log_file="${output_prefix}.histXX.withdel${pDel}.log"

# Keep terminal output unchanged, but normalize carriage returns in the log
# so progress updates do not show up as long ^M-filled lines.
exec > >(tee >(perl -pe 's/\r/\n/g' > "$log_file")) 2>&1

echo "########################################"
echo "COALESCENT SIMULATION FOR K-MER SPRECTRA WITH TERMINAL DELETION"
echo "Started: $(date '+%Y-%m-%d %H:%M:%S %Z')"

# make a temp dir
td=$(mktemp -d "${tmp_base%/}/csks.XXXXXX")
echo "Temporary directory: $td"

# run msprime (this produces the genomes)
python code/makeGenomesDel.py "$gs" "$ploy" "$theta" "$td" "$pDel"

# make reads from genomes
python code/makeReads.py "$td"

# run kmc, loop over k-mer lengths
for kk in 21 24 27 30 33; do
    hist_file="${output_prefix}.hist${kk}"
    hist_no0_file="${output_prefix}.hist${kk}.withdel${pDel}.no0"

    echo "########################################"
    echo "k = $kk"
    echo "Started k-mer analysis: $(date '+%Y-%m-%d %H:%M:%S %Z')"
    echo "########################################"
    echo "KMC DB"
    kmc -k"$kk" -t3 -cs500000000 -fm "$td/reads.fa" "$td/db$kk" "$td"
    echo "DB done."
    echo ""
    echo "########################################"
    echo "MAKING SPECTRUM"
    # run kmc_tools to make spectrum
    kmc_tools transform "$td/db$kk" histogram "$hist_file" -cx5000
    echo "removing zero lines"
    {
        printf "#gs=%s,ploy=%s,theta=%s,k=%s,pDel=%s\n" "$gs" "$ploy" "$theta" "$kk" "$pDel"
        awk '{ if( $2 != 0 ){ print $0 } }' "$hist_file"
    } > "$hist_no0_file" && rm "$hist_file"
    echo "Histogram written to $hist_no0_file"
    echo "Finished k-mer analysis: $(date '+%Y-%m-%d %H:%M:%S %Z')"
    echo ""
done

echo "########################################"
echo "Removing temp files..."
if [ "$keep_tmp" -eq 1 ]; then
    echo "Keeping temp directory at $td"
else
    rm -rf "$td"
fi
echo "Histograms written to ${output_prefix}.histXX.withdel${pDel}.no0"
echo "Log written to $log_file"
echo "Finished: $(date '+%Y-%m-%d %H:%M:%S %Z')"
echo "PIPELINE DONE."
