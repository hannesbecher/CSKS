# generate k-mer spectrum according to genome and population parameters

# activate conda environment (run 00setUpEnv.sh to set up)
conda activate CSKS

# make a temp dir 
td=$(mktemp -d)


# simulation parameters
gs=1000000 # (haploid) genome size
ploy=2 # ploidy level
theta=0.005 # populatoin-scaled mutation rate (= heterozygosity in a random-mating population)
#gc=0.5 #  genome GC-content (ignored for now)
#div=0.05 # genome divergence (allotetraploids only)

# run msprime (this produces the genomes)
# python code/makeGenomes.py $gs $ploy $theta

# make reads from genomes

# run kmc

# run kmc_tools to make spectrum


