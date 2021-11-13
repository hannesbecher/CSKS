# run msprime to generate genome sequences
# USAGE: python makeGenomes.py gs ploy theta
import msprime
import sys
import numpy as np
import random


rng = np.random.default_rng()

gc=0.5
gs, ploy, theta, tdir = sys.argv[1:5]

# tdir = "tmp/"
print("Simulating data for haploid GS %d, ploidy %d, and theta %1.3f." % (gs, ploy, theta))


chrs=10 # number of replicate runs corresponding to "chromosomes" of equal size
gs, ploy, theta = 10000, 2, 0.005
ts = msprime.sim_ancestry(
        samples=ploy,
        ploidy=1,
        population_size=1,
        discrete_genome=True,# default setting now?
        recombination_rate=1e-4, # possibly adjust
        sequence_length=gs,
        num_replicates=chrs,
        model=msprime.BinaryMutationModel) 
#        random_seed=123456)

mutated_ts = [msprime.sim_mutations(i, rate=theta) for i in ts]
#        random_seed=123456)

# array of genotypes (concatenating all chromosomes)
gts = np.zeros((chrs*gs, ploy), dtype="int8") # We only have ones and zeros. Could thus use bool array. But bool behaves different to int when used as index.


# A list of lists, one for each "chromosome"
# Each chromsome list contains tuples of variant positions and genotypes
tupListList = [[(i.site.position, i.genotypes) for i in j.variants()] for j in mutated_ts]

# A function to update the genotype array with the variant information
def updateGts(gtArr, tup, chrNo):
    gtArr[int(tup[0]) + chrNo * gs,:] = tup[1]

# Update the genotype array
for i in range(chrs):
    for j in tupListList[i]:
        updateGts(gts, j, i)

# 1D array of nucleotides:
refArr = np.random.choice(["A","C","G","T"], chrNo*gs, replace=True, p=[(1-gc)/2,(1-gc)/2, (gc)/2, (gc)/2 ])

def pickAlt(ref, gc=gc):
    a=rng.choice(["A","T","C","G"],1,replace=False,p=[(1-gc)/2,(1-gc)/2, (gc)/2, (gc)/2 ])[0]
    if a == ref:
        a = pickAlt(ref, gc=gc)
    return a

#print("Genrating random alleles...", end="")
#alleles = np.vstack(pickAlleles() for _ in range(chrs*gs))
print("Genrating alt alleles...", end="")
altDict = {i:pickAlt(refArr[i]) for i in np.where(np.sum(gts, 1)>0)}
print("Done.")

genomes = [alleles[range(chrs*gs),gts[:,i]] for i in range(ploy)]

outhandle = open(tdir + "genomes.fa", "w")
outhandle.write(">genomes\n")
for i in genomes:
    outhandle.write("".join(i))

outhandle.close()


# # a dictionary of variant site positions and corresponding alleles
# alleles = {i: pickAlleles() for i in np.where(np.sum(gts, 1)>0)[0]}





#len(np.where(gts[:,0] > 0)[0])


