# activate 'conda' environment (run 00setUpEnv.sh to set up)
#micromamba activate CSKS

usage() {
    cat <<EOF
Usage: $0 [-g genome_size] [-p ploidy] [-t theta] [-k kmer_length] [-d temp_base] [--keep] [-h]

Run the CSKS simulation pipeline and write a k-mer histogram.

Options:
  -g, --gs       Chromosome size in nt (default: 1000000)
  -p, --ploy     Ploidy level (default: 2)
  -t, --theta    Population-scaled mutation rate (default: 0.02)
  -k, --k        K-mer length for KMC analysis (default: 21)
  -d, --tmpdir   Base path for the temporary directory (default: /tmp)
  --keep         Keep the temporary directory for debugging
  -h, --help     Show this help message and exit
EOF
}

# simulation parameters
gs=1000000 # Size in nt of each chromosome. There are 10 chromosomes, so the overall haploid genome size is 'gs'*10!
ploy=2 # ploidy level
theta=0.02 # population-scaled mutation rate (= heterozygosity in a random-mating population)
k=21 # k-mer length
tmp_base=/tmp # base path for temporary directory
keep_tmp=0 # keep temporary directory for debugging
#gc=0.5 #  genome GC-content (ignored for now)
#div=0.05 # genome divergence (allotetraploids only)

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
        -k|--k)
            if [ -z "$2" ]; then
                echo "Missing value for $1" >&2
                usage >&2
                exit 1
            fi
            k="$2"
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

hist_file="G${gs}P${ploy}T${theta}.hist"
hist_no0_file="${hist_file}.no0"
log_file="G${gs}P${ploy}T${theta}.log"

# Keep terminal output unchanged, but normalize carriage returns in the log
# so progress updates do not show up as long ^M-filled lines.
exec > >(tee >(perl -pe 's/\r/\n/g' > "$log_file")) 2>&1

# generate k-mer spectrum according to genome and population parameters
echo "########################################"
echo "COALESCENT SIMULATION FOR K-MER SPRECTRA"
echo "Started: $(date '+%Y-%m-%d %H:%M:%S %Z')"

# make a temp dir
td=$(mktemp -d "${tmp_base%/}/csks.XXXXXX")

# run msprime (this produces the genomes)
python code/makeGenomes.py "$gs" "$ploy" "$theta" "$td"

# make reads from genomes
python code/makeReads.py "$td"

# run kmc
echo "########################################"
echo "KMC DB"
echo "Started: $(date '+%Y-%m-%d %H:%M:%S %Z')"
kmc -k"$k" -t3 -cs500000000 -fm $td/reads.fa $td/db21 $td
echo "DB done."
echo ""
echo "########################################"
echo "MAKING SPECTRUM"
echo "Started: $(date '+%Y-%m-%d %H:%M:%S %Z')"
# run kmc_tools to make spectrum
kmc_tools transform  $td/db21 histogram  "$hist_file" -cx5000
echo "removing zero lines"
{
    printf "#gs=%E,ploy=%s,theta=%s,k=%s\n" "$gs" "$ploy" "$theta" "$k"
    awk '{ if( $2 != 0 ){ print $0 } }' "$hist_file"
} > "$hist_no0_file" && rm "$hist_file"
echo "Done."
echo ""

echo "########################################"
echo "Removing temp files..."
echo "Started: $(date '+%Y-%m-%d %H:%M:%S %Z')"
if [ "$keep_tmp" -eq 1 ]; then
    echo "Keeping temp directory at $td"
else
    rm -rf "$td"
fi
echo "Histogram written to $hist_no0_file"
echo "Log written to $log_file"
echo "PIPELINE DONE."
echo "At: $(date '+%Y-%m-%d %H:%M:%S %Z')"