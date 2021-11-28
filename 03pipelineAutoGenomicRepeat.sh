# generate k-mer spectrum according to genome and population parameters
echo "########################################"
echo "COALESCENT SIMULATION FOR K-MER SPRECTRA"

# activate conda environment (run 00setUpEnv.sh to set up)
conda activate CSKS

# make a temp dir 
td=$(mktemp -d)

suf=$1
# simulation parameters
gs=1000000 # Size in nt of each chromosome. There are 10 chromosomes, so the overall haploid genome size is 'gs'*10!
#ploy=4 # ploidy level
ploy=2 # ploidy level
theta=0.005 # population-scaled mutation rate (= heterozygosity in a random-mating population)
#theta=0.01 # population-scaled mutation rate (= heterozygosity in a random-mating population)
#gc=0.5 #  genome GC-content (ignored for now)
#div=0.05 # genome divergence (allotetraploids only)

# run msprime (this produces the genomes)
#python code/makeGenomes.py $gs $ploy $theta $td
python code/makeGenomesRateMapGenomicRepeat.py $gs $ploy $theta $td
#python makeGenomes.py $gs 2 2 0.01 20 $td

# make reads from genomes
python code/makeReads.py $td

# run kmc
echo "########################################"
echo "KMC DB"
kmc -k21 -t3 -cs500000000 -fm $td/reads.fa $td/db21 $td
echo "DB done."
echo ""
echo "########################################"
echo "MAKING SPECTRUM"
# run kmc_tools to make spectrum
kmc_tools transform  $td/db21 histogram  G$gs'P'$ploy'T'$theta$suf'GenRep.hist' -cx5000
echo "removing zero lines of multiplicity > 500"
awk '{ if( ($1 < 501) || ($2 != 0 )){ print $0 } }' G$gs'P'$ploy'T'$theta$suf'GenRep.hist' > G$gs'P'$ploy'T'$theta$suf'GenRep.hist'.no0 && rm G$gs'P'$ploy'T'$theta$suf'GenRep.hist'
echo "Done."
echo ""

echo "########################################"
echo "Removing temp files..."
rm -rf $td # comment out for debug!
echo "Histogram written to G"$gs"P"$ploy"T"$theta$suf"GenRep.hist.no0"
echo "PIPELINE DONE."
echo "########################################\n"
