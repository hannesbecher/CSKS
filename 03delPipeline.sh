# generate k-mer spectrum according to genome and population parameters
echo "########################################"
echo "COALESCENT SIMULATION FOR K-MER SPRECTRA"

# activate conda environment (run 00setUpEnv.sh to set up)
#conda activate CSKS_nomap

# make a temp dir 
td=$(mktemp -d)


# simulation parameters
gs=100000 # Size in nt of each chromosome. There are 10 chromosomes, so the overall haploid genome size is 'gs'*10!
# ploidy level is 2
theta=0.02 # population-scaled mutation rate (= heterozygosity in a random-mating population)
gc=0.5 #  genome GC-content (ignored for now)
pDel=0.1 # prop deleted

# run genome script
python code/makeDiploidNoCoalDel.py $gs $theta $gc $td $pDel

# make reads from genomes
python code/makeReads.py $td

# run kmc
echo "########################################"
echo "KMC DB"
kmc -k21 -t3 -cs500000000 -fm $td/reads.fa $td/db21 $td
kmc -k24 -t3 -cs500000000 -fm $td/reads.fa $td/db24 $td
kmc -k27 -t3 -cs500000000 -fm $td/reads.fa $td/db27 $td
kmc -k30 -t3 -cs500000000 -fm $td/reads.fa $td/db30 $td
kmc -k33 -t3 -cs500000000 -fm $td/reads.fa $td/db33 $td
echo "DB done."
echo ""
echo "########################################"
echo "MAKING SPECTRUM"
# run kmc_tools to make spectrum
kmc_tools transform  $td/db21 histogram  G$gs'P'$ploy'T'$theta.hist21 -cx5000
kmc_tools transform  $td/db24 histogram  G$gs'P'$ploy'T'$theta.hist24 -cx5000
kmc_tools transform  $td/db27 histogram  G$gs'P'$ploy'T'$theta.hist27 -cx5000
kmc_tools transform  $td/db30 histogram  G$gs'P'$ploy'T'$theta.hist30 -cx5000
kmc_tools transform  $td/db33 histogram  G$gs'P'$ploy'T'$theta.hist33 -cx5000
echo "removing zero lines"
awk '{ if( $2 != 0 ){ print $0 } }' G$gs'P'$ploy'T'$theta.hist21 > G$gs'P'$ploy'T'$theta.hist21.no0 && rm G$gs'P'$ploy'T'$theta.hist21
awk '{ if( $2 != 0 ){ print $0 } }' G$gs'P'$ploy'T'$theta.hist24 > G$gs'P'$ploy'T'$theta.hist24.no0 && rm G$gs'P'$ploy'T'$theta.hist24
awk '{ if( $2 != 0 ){ print $0 } }' G$gs'P'$ploy'T'$theta.hist27 > G$gs'P'$ploy'T'$theta.hist27.no0 && rm G$gs'P'$ploy'T'$theta.hist27
awk '{ if( $2 != 0 ){ print $0 } }' G$gs'P'$ploy'T'$theta.hist30 > G$gs'P'$ploy'T'$theta.hist30.no0 && rm G$gs'P'$ploy'T'$theta.hist30
awk '{ if( $2 != 0 ){ print $0 } }' G$gs'P'$ploy'T'$theta.hist33 > G$gs'P'$ploy'T'$theta.hist33.no0 && rm G$gs'P'$ploy'T'$theta.hist33
echo "Done."
echo ""

echo "########################################"
echo "Removing temp files..."
rm -rf $td # comment out for debug!
echo "Histogram written to G"$gs"P"$ploy"T"$theta".histXX.no0"
echo "PIPELINE DONE."

