# generate k-mer spectrum according to genome and population parameters
echo "########################################"
echo "COALESCENT SIMULATION FOR K-MER SPRECTRA"

# activate conda environment (run 00setUpEnv.sh to set up)
#conda activate CSKS

# make a temp dir 
td=$(mktemp -d)


# simulation parameters
gs=1000000 # Size in nt of each chromosome. There are 10 chromosomes, so the overall haploid genome size is 'gs'*10!
ploy=2 # ploidy level
theta=0.002 # population-scaled mutation rate (= heterozygosity in a random-mating population)
pDel=0.001 # proportion deleted
#gc=0.5 #  genome GC-content (ignored for now)
#div=0.05 # genome divergence (allotetraploids only)

# run msprime (this produces the genomes)
python code/makeGenomesDel.py $gs $ploy $theta $td $pDel # gc is not passed ATM



# make reads from genomes
python code/makeReads.py $td



# run kmc, loop over k-mer lengths
for kk in 21 24 27 30 33; do
	echo "########################################"
	echo "k = $kk."
	echo "########################################"
	echo "KMC DB"
	kmc -k$kk -t3 -cs500000000 -fm $td/reads.fa $td/db$kk $td
	echo "DB done."
	echo ""
	echo "########################################"
	echo "MAKING SPECTRUM"
	# run kmc_tools to make spectrum
	kmc_tools transform  $td/db$kk histogram  G$gs'P'$ploy'T'$theta.hist$kk -cx5000
	echo "removing zero lines"
	awk '{ if( $2 != 0 ){ print $0 } }' G$gs'P'$ploy'T'$theta.hist$kk > G$gs'P'$ploy'T'$theta.hist$kk.no0 && rm G$gs'P'$ploy'T'$theta.hist$kk
	echo "Done."
	echo ""
	
done

echo "########################################"
echo "Removing temp files $td..."
rm -rf $td # comment out for debug!
echo "Histograms written to G"$gs"P"$ploy"T"$theta".histXX.no0"
echo "PIPELINE DONE."
