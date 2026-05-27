# CSKS
Coalescent simulations of k-mer spectra for testing k-mer genome profiling tools

## Usage
1. Set up (conda) environment, see example commands in `00setUpEnv`
2. Activate the newly created environment, e.g. `conda activate CSKS` 
3. Run pipeline, e.g.: `bash 05delPipeCoalComplex.sh  -t 0.01 -g 10000001 -p 2 --pdel 0.0 --error-profile errProfile`
4. Check helpstring `bash 05delPipeCoalComplex.sh -h` for more options and details on parameters.
4. Adjust parameters and re-run if desired.


## The motivation
Simulating genetic data is not trivial. In the case of one diploid (i.e. two genomes), one can just throw in heterozygous sites according to the desired heterozygosity level (θ). If more than two genomes are simulated, the task becomes more complicated. This is because related genomes have a common ancestor and thus share mutations. A summary of the number of times different mutations are shared between genomes in a sample is called the [site-frequency spectrum](https://en.wikipedia.org/wiki/Allele_frequency_spectrum).

## Implementation
There are at the moment four steps that take place in a temporary directory. Only the final result, a k-mer spectrum, is stored in the working directory. The steps are:
1. The simulation of genetic variants and storage of simulated genomes
2. Generation of sequencing reads from the simulated genomes
3. Generation of a or multiple KMC3 database(s)
4. Generation of a k-mer spectrum or multiple spectra from the KMC3 database(s)

Below, these steps are descibed in more detail for the standard entry point, which currently the script `05delPipeCoalComplex.sh`.

### Coalescent simulation
This is done by the script `makeGenomes.py`, which runs the coalescent simulator [`msprime`](https://tskit.dev/msprime/docs/stable/quickstart.html) to generate genetic variant data in a sample of `n` genomes (where `n` is the ploidy level specified). `msprime` generates information only for variant sites and so the invariant rest of the genome is made up separately.

The genomes are then written (concatenated) to a single-record FASTA file, `genomes.fa`, with a specific header format containing the genome size information: `>genomes_p[ploidy]_s[haploid GS]_t[total GS]`.

This script can be modified to allow for allopolyploidy, genomic repeats, unequal sub-genomes, etc.

### Generation of sequencing reads
This is done by `makeReadsLowComplex.py`. This script slices reads out of the FASTA reference file. It can add sequencing errors based on a constant error rate (`--error-prob`) or based on a profile file listing one error rate for each position in the read (`--error-profile`). Only constand read length are supported ATM. The reads are written to `reads.fa`.

### K-mer data base
KMC3 is run on `reads.fa` with k-mer lengths ranging from 21 to 51 in steps of 3.

### Spectra
11 k-mer spectra are generated from the KMC DBs and the temporary files are deleted by default (can be kept with `--keep`).

