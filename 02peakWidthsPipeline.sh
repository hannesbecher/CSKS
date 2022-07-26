# generate k-mer spectrum according to genome and population parameters
echo "########################################"
echo "COALESCENT SIMULATION FOR K-MER SPRECTRA"

# activate conda environment (run 00setUpEnv.sh to set up)
#conda activate CSKS

# make a temp dir 
#td=$(mktemp -d)
td=/tmp/tmp.PwXQCJ1tQi

echo "TMP DIR $td"

# simulation parameters
gs=1000000 # Size in nt of each chromosome. There are 10 chromosomes, so the overall haploid genome size is 'gs'*10!
#chrs=1
#ploy=2 # ploidy level
theta=0.001 # population-scaled mutation rate (= heterozygosity in a random-mating population)
gc=0.4 #  genome GC-content (ignored for now)
#div=0.05 # genome divergence (allotetraploids only)
rl=150 # read length
dep=100



# run msprime (this produces the genomes)
python code/makeDiploidNoCoal.py $gs $theta $gc $td

# make reads from genomes
#python code/makeReads.py $td
wgsim -r 0 -R 0 -e 0.001 -d 550 -s 50 -1 $rl -2 $rl -N $((gs/2/$rl*$dep)) $td/genomes.fa $td/fw.fq $td/rw.fq > /dev/null


# run kmc
echo "########################################"
echo "KMC DB"
printf "$td/fw.fq\n$td/rw.fq" > $td/readFiles
kmc -k21 -t3 -cs500000000 -fq @$td/readFiles $td/db21 $td
echo "DB done."
echo ""

echo "########################################"
echo "MAKING SPECTRUM"
# run kmc_tools to make spectrum
kmc_tools transform  $td/db21 histogram  G$gs'P2T'$theta.hist -cx5000
echo "removing zero lines"
awk '{ if( $2 != 0 ){ print $0 } }' G$gs'P2T'$theta.hist > G$gs'P2T'$theta.hist.no0 && rm G$gs'P2T'$theta.hist
echo "Done."
echo ""

echo "########################################"
echo "INDEXING REFERENCE"
bwa-mem2 index $td/reference.fa

echo "########################################"
echo "MAPPING"
bwa-mem2 mem -t 6 $td/reference.fa $td/fw.fq $td/rw.fq | samtools view -b -o $td/ref.bam

echo "########################################"
echo "SORTING BAM"
samtools sort -o $td/refS.bam $td/ref.bam && rm $td/ref.bam

echo "########################################"
echo "MARKING DUPLICATES"
picard MarkDuplicates I=$td/refS.bam O=$td/refSD.bam M=$td/refSD.dup-metrics.txt && rm $td/refS.bam

echo "########################################"
echo "EXTRACTING DUPLICATES"
samtools view -f 1024 -b -o $td/refSDonly.bam $td/refSD.bam
samtools index $td/refSDonly.bam
samtools index $td/refSD.bam

echo "########################################"
echo "COMPUTING MAPPING DEPTHS"
mosdepth -b 1000 $td/refSD $td/refSD.bam 
mosdepth -b 1000 $td/refSDonly $td/refSDonly.bam 



#echo "########################################"
#echo "Removing temp files..."
#rm -rf $td # comment out for debug!
#echo "Histogram written to G"$gs"P"$ploy"T"$theta".hist.no0"
#echo "PIPELINE DONE."

