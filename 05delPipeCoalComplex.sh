# activate conda environment (run 00setUpEnv.sh to set up)
#conda activate CSKS

script_dir=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

usage() {
    cat <<EOF
Usage: $0 [-g genome_size] [-p ploidy] [-t theta] [--pdel proportion_deleted] [--cplx library_complexity] [--read-length read_length] [--depth sequencing_depth] [--error-prob probability | --error-profile file] [--pref output_prefix] [-d temp_base] [--keep] [-h]

Run the coalescent deletion pipeline and write k-mer histograms for multiple k values.

Options:
  -g, --gs       Chromosome size in nt (default: 10000000)
  -p, --ploy     Ploidy level (default: 2)
  -t, --theta    Population-scaled mutation rate (default: 0.1)
  --pdel         Proportion deleted (default: 0.01)
  --cplx         Library complexity (0.0 to 1.0, default: 1.0)
  --read-length  Read length in nt (default: 150)
  --depth        Sequencing coverage per haploid genome (default: 50)
  --error-prob   Fixed per-nucleotide sequencing error probability (default: 0.0)
  --error-profile
                 Text file with one error probability per read position
  --pref         Output file prefix (default: G[gs]P[ploidy]T[theta]C[cplx])
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
cplx=1.0 # library complexity
read_length=150 # read length
depth=50 # sequencing coverage per haploid genome
error_prob=0.0 # fixed per-nucleotide sequencing error probability
error_profile= # file with per-position sequencing error probabilities
error_mode=none # one of: none, fixed, profile
output_prefix= # output file prefix; default is set after parsing
tmp_base=/tmp # base path for temporary directory
keep_tmp=0 # keep temporary directory for debugging
kmer_start=21 # first k-mer length
kmer_step=3 # k-mer length step
kmer_end=51 # last k-mer length
kmc_threads=3 # threads used by KMC
kmc_counter_max=500000000 # maximum KMC counter value
hist_cutoff_max=5000 # maximum count included in KMC histogram output

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
        --cplx)
            if [ -z "$2" ]; then
                echo "Missing value for $1" >&2
                usage >&2
                exit 1
            fi
            cplx="$2"
            # Validate cplx is a float between 0 and 1
            if ! awk -v val="$cplx" 'BEGIN { if (val < 0 || val > 1 || val !~ /^[0-9]*\.?[0-9]+$/) exit 1 }'; then
                echo "cplx must be a float between 0.0 and 1.0" >&2
                usage >&2
                exit 1
            fi
            shift 2
            ;;
        --read-length|--rlen)
            if [ -z "$2" ]; then
                echo "Missing value for $1" >&2
                usage >&2
                exit 1
            fi
            read_length="$2"
            if ! awk -v val="$read_length" 'BEGIN { if (val !~ /^[0-9]+$/ || val <= 0) exit 1 }'; then
                echo "read length must be a positive integer" >&2
                usage >&2
                exit 1
            fi
            shift 2
            ;;
        --depth|--coverage)
            if [ -z "$2" ]; then
                echo "Missing value for $1" >&2
                usage >&2
                exit 1
            fi
            depth="$2"
            if ! awk -v val="$depth" 'BEGIN { if (val !~ /^[0-9]+$/ || val <= 0) exit 1 }'; then
                echo "depth must be a positive integer" >&2
                usage >&2
                exit 1
            fi
            shift 2
            ;;
        --error-prob)
            if [ -z "$2" ]; then
                echo "Missing value for $1" >&2
                usage >&2
                exit 1
            fi
            if [ "$error_mode" != "none" ]; then
                echo "Use either --error-prob or --error-profile, not both" >&2
                usage >&2
                exit 1
            fi
            error_prob="$2"
            if ! awk -v val="$error_prob" 'BEGIN { if (val < 0 || val > 1 || val !~ /^[0-9]*\.?[0-9]+$/) exit 1 }'; then
                echo "error probability must be a float between 0.0 and 1.0" >&2
                usage >&2
                exit 1
            fi
            error_mode=fixed
            shift 2
            ;;
        --error-profile)
            if [ -z "$2" ]; then
                echo "Missing value for $1" >&2
                usage >&2
                exit 1
            fi
            if [ "$error_mode" != "none" ]; then
                echo "Use either --error-prob or --error-profile, not both" >&2
                usage >&2
                exit 1
            fi
            error_profile="$2"
            if [ ! -f "$error_profile" ]; then
                echo "error profile file does not exist: $error_profile" >&2
                usage >&2
                exit 1
            fi
            error_mode=profile
            shift 2
            ;;
        --pref)
            if [ -z "$2" ]; then
                echo "Missing value for $1" >&2
                usage >&2
                exit 1
            fi
            output_prefix="$2"
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

if [ -z "$output_prefix" ]; then
    output_prefix="G${gs}P${ploy}T${theta}C${cplx}"
fi
log_file="${output_prefix}.histXX.withdel${pDel}.log"

# Keep terminal output unchanged, but normalize carriage returns in the log
# so progress updates do not show up as long ^M-filled lines.
exec > >(tee >(perl -pe 's/\r/\n/g' > "$log_file")) 2>&1

echo "########################################"
echo "COALESCENT SIMULATION FOR K-MER SPRECTRA WITH TERMINAL DELETION"
echo "Started: $(date '+%Y-%m-%d %H:%M:%S %Z')"
echo "Input parameters:"
echo "  Script directory: $script_dir"
echo "  Chromosome size: $gs"
echo "  Ploidy: $ploy"
echo "  Theta: $theta"
echo "  Proportion deleted: $pDel"
echo "  Library complexity: $cplx"
echo "  Read length: $read_length"
echo "  Sequencing depth: $depth"
echo "  Error mode: $error_mode"
if [ -n "$error_profile" ]; then
    echo "  Error profile: $error_profile"
else
    echo "  Error probability: $error_prob"
fi
echo "  Temporary directory base: $tmp_base"
echo "  Keep temporary directory: $keep_tmp"
echo "  Output file prefix: $output_prefix"
echo "  Log file: $log_file"
echo "  K-mer lengths: ${kmer_start}..${kmer_end} step ${kmer_step}"
echo "  KMC threads: $kmc_threads"
echo "  KMC counter max: $kmc_counter_max"
echo "  Histogram cutoff max: $hist_cutoff_max"

# make a temp dir
td=$(mktemp -d "${tmp_base%/}/csks.XXXXXX")
echo "Temporary directory: $td"

# run msprime (this produces the genomes)
python "$script_dir/code/makeGenomesDel.py" "$gs" "$ploy" "$theta" "$td" "$pDel"

# make reads from genomes
read_error_args=(--error-prob "$error_prob")
if [ -n "$error_profile" ]; then
    read_error_args=(--error-profile "$error_profile")
fi
python "$script_dir/code/makeReadsLowComplex.py" "$td" "$cplx" --read-length "$read_length" --depth "$depth" "${read_error_args[@]}"

# run kmc, loop over k-mer lengths
for kk in `seq "$kmer_start" "$kmer_step" "$kmer_end"`; do
    hist_file="${output_prefix}.hist${kk}"
    hist_no0_file="${output_prefix}.hist${kk}.withdel${pDel}.no0"

    echo "########################################"
    echo "k = $kk"
    echo "Started k-mer analysis: $(date '+%Y-%m-%d %H:%M:%S %Z')"
    echo "########################################"
    echo "KMC DB"
    kmc -k"$kk" -t"$kmc_threads" -cs"$kmc_counter_max" -fm "$td/reads.fa" "$td/db$kk" "$td"
    echo "DB done."
    echo ""
    echo "########################################"
    echo "MAKING SPECTRUM"
    # run kmc_tools to make spectrum
    kmc_tools transform "$td/db$kk" histogram "$hist_file" -cx"$hist_cutoff_max"
    echo "removing zero lines"
    {
        printf "#gs=%s,ploy=%s,theta=%s,k=%s,pDel=%s,cplx=%s\n" "$gs" "$ploy" "$theta" "$kk" "$pDel" "$cplx"
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
