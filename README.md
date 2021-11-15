# CSKS
Coalescent simulations the generate k-mer spectra

## Usage
1. Set up conda environment: `bash 00setUpEnv.sh`
2. Run pipeline: `bash -i 01pipeline.sh`
3. Adjust paprameters in `01pipeline.sh` and re-run if desired.


## The motivation
Simulating genetic data is not trivial. In the case of one diploid (i.e. two genomes), one can just throw in heterozygous sites according to the desired heterozygosity level (θ). If more than two genomes are simulated, the task becomes more complicated. This is because related genomes have a common ancestor and thus share mutations. A summary of the number of times different mutations are shared between genomes in a sample is called the [site-frequency spectrum](https://en.wikipedia.org/wiki/Allele_frequency_spectrum).

## Implementation
There are at the moment four steps that take place in a temporary directory. Only the final result, a k-mer spectrum, is stored in the working directory. The steps are:
1. The simulation of genetic variants and storage of simulated genomes
2. Generation of sequencing reads from the simulated genomes
3. Generation of a KMC3 database
4. Generation of a k-mer spectrum

Below, these steps are descibed in more detail.

### Coalescent simulation
This is done by the script `makeGenomes.py`, which runs the coalescent simulator [`msprime`](https://tskit.dev/msprime/docs/stable/quickstart.html) to generate genetic variant data in a sample of `n` genomes (where `n` is the ploidy level specified). `msprime` generates information only for variant sites and so the invariant rest of the genome is made up separately.

The genomes are then written (concatenated) to a single-record FASTA file, `genomes.fa`, with a specific header format containing the genome size information: `>genomes_p[ploidy]_s[haploid GS]_t[total GS]`.

This script can be modified to allow for allopolyploidy, genomic repeats, unequal sub-genomes, etc.

### Generation of sequencing reads
This is done by `makeReads.py`, a very simple script that may well be replaced in the future by a more sophisticated sequnecing read simulator (there are several available for paired reads with error models).

For now, this generates error-free single-end reads in FASTA format which are written to `reads.fa`. To retrieve the genome size, the script looks at the header line of `genomes.fa`. It then generates reads of 150nt at 30x sequencing coverage. The read positions are sampled at random from a uniform distribution on the interval \[1, genome size-150\]. This should lead to narrow 'Poisson' peaks (bias parameter 0).

### K-mer data base
KMC3 is run on `reads.fa` with the k-mer length set to 21.

### Spectra
A k-mer spectrum is generated from the KMC DB and the temporary files are deleted.

