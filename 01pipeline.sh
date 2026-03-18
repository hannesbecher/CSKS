# generate k-mer spectrum according to genome and population parameters
echo "########################################"
echo "COALESCENT SIMULATION FOR K-MER SPRECTRA"

# activate 'conda' environment (run 00setUpEnv.sh to set up)
#micromamba activate CSKS

usage() {
    cat <<EOF
Usage: $0 [-g genome_size] [-p ploidy] [-t theta] [-k kmer_length] [-h]

Run the CSKS simulation pipeline and write a k-mer histogram.

Options:
  -g, --gs       Chromosome size in nt (default: 1000000)
  -p, --ploy     Ploidy level (default: 2)
  -t, --theta    Population-scaled mutation rate (default: 0.02)
  -k, --k        K-mer length for KMC analysis (default: 21)
  -h, --help     Show this help message and exit
EOF
}

# make a temp dir 
td=$(mktemp -d)


# simulation parameters
gs=1000000 # Size in nt of each chromosome. There are 10 chromosomes, so the overall haploid genome size is 'gs'*10!
ploy=2 # ploidy level
theta=0.02 # population-scaled mutation rate (= heterozygosity in a random-mating population)
k=21 # k-mer length
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

# run msprime (this produces the genomes)
python code/makeGenomes.py "$gs" "$ploy" "$theta" "$td"

# make reads from genomes
python code/makeReads.py "$td"

# run kmc
echo "########################################"
echo "KMC DB"
kmc -k"$k" -t3 -cs500000000 -fm $td/reads.fa $td/db21 $td
echo "DB done."
echo ""
echo "########################################"
echo "MAKING SPECTRUM"
# run kmc_tools to make spectrum
kmc_tools transform  $td/db21 histogram  G$gs'P'$ploy'T'$theta.hist -cx5000
echo "removing zero lines"
awk '{ if( $2 != 0 ){ print $0 } }' G$gs'P'$ploy'T'$theta.hist > G$gs'P'$ploy'T'$theta.hist.no0 && rm G$gs'P'$ploy'T'$theta.hist
echo "Done."
echo ""

echo "########################################"
echo "Removing temp files..."
rm -rf "$td" # comment out for debug!
echo "Histogram written to G"$gs"P"$ploy"T"$theta".hist.no0"
echo "PIPELINE DONE."
