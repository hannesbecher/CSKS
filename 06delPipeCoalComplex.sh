# activate conda environment (run 00setUpEnv.sh to set up)
#conda activate CSKS

script_dir=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

usage() {
    cat <<EOF
Usage: $0 [-g genome_size] [-p ploidy] [-t theta] [--pdel proportion_deleted] [--cplx library_complexity] [--read-length read_length] [--depth sequencing_depth] [--error-prob probability | --error-profile file] [--pcr-dist none|geometric|lognormal] [--pcr-param P] [--pref output_prefix] [-d temp_base] [--keep] [-h]

Run the coalescent deletion pipeline with an imperfect-PCR read simulation
step and write k-mer histograms for multiple k values.

This is a copy of 05delPipeCoalComplex.sh that calls
code/makeReadsImperfectPCR.py instead of code/makeReadsLowComplex.py. All
options are identical to 05delPipeCoalComplex.sh, with two new ones:

  --pcr-dist     Distribution of per-molecule PCR amplification factors:
                 none (default; reproduces the original 05... behaviour
                 exactly), geometric, or lognormal
  --pcr-param    Parameter of the amplification-factor distribution.
                 geometric: success probability p in (0, 1] (mean copies
                 = 1/p; p=1 means no amplification noise).
                 lognormal: log-scale standard deviation sigma >= 0
                 (sigma=0 means no amplification noise).
                 Required unless --pcr-dist none.

Other options:
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
  --pref         Output file prefix (default: G[gs]P[ploidy]T[theta]C[cplx]PCR[dist][param])
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
pcr_dist=none # one of: none, geometric, lognormal
pcr_param= # parameter of the amplification-factor distribution
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
        --pcr-dist)
            if [ -z "$2" ]; then
                echo "Missing value for $1" >&2
                usage >&2
                exit 1
            fi
            pcr_dist="$2"
            case "$pcr_dist" in
                none|geometric|lognormal) ;;
                *)
                    echo "--pcr-dist must be one of: none, geometric, lognormal" >&2
                    usage >&2
                    exit 1
                    ;;
            esac
            shift 2
            ;;
        --pcr-param)
            if [ -z "$2" ]; then
                echo "Missing value for $1" >&2
                usage >&2
                exit 1
            fi
            pcr_param="$2"
            if ! awk -v val="$pcr_param" 'BEGIN { if (val !~ /^[0-9]*\.?[0-9]+$/) exit 1 }'; then
                echo "--pcr-param must be a non-negative number" >&2
                usage >&2
                exit 1
            fi
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

if [ "$pcr_dist" != "none" ] && [ -z "$pcr_param" ]; then
    echo "--pcr-param is required when --pcr-dist is $pcr_dist" >&2
    usage >&2
    exit 1
fi
if [ "$pcr_dist" = "none" ] && [ -n "$pcr_param" ]; then
    echo "--pcr-param must not be set when --pcr-dist is none" >&2
    usage >&2
    exit 1
fi

if [ -z "$output_prefix" ]; then
    if [ "$pcr_dist" = "none" ]; then
        output_prefix="G${gs}P${ploy}T${theta}C${cplx}PCRnone"
    else
        output_prefix="G${gs}P${ploy}T${theta}C${cplx}PCR${pcr_dist}${pcr_param}"
    fi
fi
log_file="${output_prefix}.histXX.withdel${pDel}.log"

# Keep terminal output unchanged, but normalize carriage returns in the log
# so progress updates do not show up as long ^M-filled lines.
exec > >(tee >(perl -pe 's/\r/\n/g' > "$log_file")) 2>&1

echo "########################################"
echo "COALESCENT SIMULATION FOR K-MER SPECTRA WITH TERMINAL DELETION AND IMPERFECT PCR"
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
echo "  PCR amplification distribution: $pcr_dist"
if [ "$pcr_dist" != "none" ]; then
    echo "  PCR amplification parameter: $pcr_param"
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

# make reads from genomes, with imperfect PCR amplification
read_error_args=(--error-prob "$error_prob")
if [ -n "$error_profile" ]; then
    read_error_args=(--error-profile "$error_profile")
fi
pcr_args=(--pcr-dist "$pcr_dist")
if [ "$pcr_dist" != "none" ]; then
    pcr_args+=(--pcr-param "$pcr_param")
fi
python "$script_dir/code/makeReadsImperfectPCR.py" "$td" "$cplx" --read-length "$read_length" --depth "$depth" "${read_error_args[@]}" "${pcr_args[@]}"

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
        printf "#gs=%s,ploy=%s,theta=%s,k=%s,pDel=%s,cplx=%s,pcr_dist=%s,pcr_param=%s\n" "$gs" "$ploy" "$theta" "$kk" "$pDel" "$cplx" "$pcr_dist" "$pcr_param"
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
