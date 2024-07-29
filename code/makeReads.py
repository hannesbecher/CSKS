# make sequencing reads from genomes file (generated using makeGenomes.py)

# USAGE: python makeReads.py tdir
import sys
import numpy as np

print("########################################\nREAD MAKING SCRIPT")

tdir = sys.argv[1]

print("Reading genomes...")
# the file has two lines: header and data
with open(tdir + "/genomes.fa", "r") as infile:
    hd, dat = infile.readlines() # there are only two lines: header and data

cov = 30 # seq coverage (per haploid genome)
rlen = 150 # read length

pl = int(hd[:-1].split("_")[1][1:])    # get ploidy from header line
hgs = int(hd[:-1].split("_")[2][1:])   # get haploid genome size from header line
tlen = int(hd[:-1].split("_")[-1][1:]) # get total genome size from header line

# compute number of reads to generate
nreads = hgs * cov * pl // rlen

# sample read locations
rstarts = np.random.choice(range(tlen-rlen), nreads, replace=True)

print("Writing reads to %s/reads.fa..." % tdir)
with open(tdir + "/reads.fa", "w") as outfile:
    c = 0
    for i in rstarts:
        outfile.write(">%d\n%s\n" % (c, dat[i:(i+rlen+1)]))
        c += 1
  
print("Reads done.\n")
