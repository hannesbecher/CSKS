# run msprime to generate genome sequences
# USAGE: python makeGenomes.py gs ploy theta
import msprime
import sys
import numpy as np
import random
import copy

print("########################################\nGENOME MAKING SCRIPT")
rng = np.random.default_rng()

gc=0.5
gs = int(sys.argv[1])   # size of one chromosome
ploy = int(sys.argv[2])
theta = float(sys.argv[3])
tdir = sys.argv[4]
pDel = float(sys.argv[5])

#tdir = "/tmp/tmp.0gVEKSPNrd"
#gs, ploy, theta = 10000, 2, 0.005
print("Simulating data for haploid GS %d, ploidy %d, and theta %1.3f..." % (gs, ploy, theta))

# one chromosome only for now
chrs=1 # number of replicate runs corresponding to "chromosomes" of equal size

ts = msprime.sim_ancestry(
        samples=ploy,
        ploidy=1,
        population_size=1,
        discrete_genome=True,# default setting now?
        recombination_rate=1e-4, # possibly adjust
        sequence_length=gs,
        num_replicates=chrs) 
#        random_seed=123456)

print("Adding mutations...")
mutated_ts = msprime.sim_mutations(i, rate=theta, model=msprime.BinaryMutationModel())
#        random_seed=123456)


asFasta = mutated_ts.as_fasta(reference_sequence=tskit.random_nucleotides(ts.sequence_length), wrap_width=0).split("\n")

# in a diploid, the 2nd genome is at list index 3


# make one deletion, at end of chr
dLen = int(pDel * gs)

asFasta[3] = asFasta[3][:-dLen]


# array of genotypes (concatenating all chromosomes)
print("Making a GT array...")
gts = np.zeros((chrs*gs, ploy), dtype=int)


# A list of lists, one for each "chromosome"
# Each chromsome list contains tuples of variant positions and genotypes
print("Extracting variant site information and population GT array...")
tupListList = [[(i.site.position, i.genotypes) for i in j.variants()] for j in mutated_ts]

# A function to update the genotype array with the variant information
def updateGts(gtArr, tup, chrNo):
    gtArr[int(tup[0]) + chrNo * gs,:] = tup[1]

# Update the genotype array
for i in range(chrs):
    for j in tupListList[i]:
        updateGts(gts, j, i)

# 1D array of nucleotides:
print("Making a random genome ref, length %d..." % (gs*chrs))
refArr = np.random.choice(["A","C","G","T"], chrs*gs, replace=True, p=[(1-gc)/2,(1-gc)/2, (gc)/2, (gc)/2 ])

def pickAlt(ref, gc=gc):
    a=rng.choice(["A","T","C","G"],1,replace=False,p=[(1-gc)/2,(1-gc)/2, (gc)/2, (gc)/2 ])[0]
    if a == ref:
        a = pickAlt(ref, gc=gc)
    return a

#print("Genrating random alleles...", end="")
#alleles = np.vstack(pickAlleles() for _ in range(chrs*gs))
print("Generating alt alleles...", end="")
altDict = {i:pickAlt(refArr[i]) for i in np.where(np.sum(gts, 1)>0)[0]}
varSites = altDict.keys()
print("Done.")

def makeGenome(k):
    if k==0: return refArr
    a = copy.deepcopy(refArr)
    varS = np.where(gts[:,k] != 0)[0]
    for i in varS:
        a[i] = altDict[i]
        
    return a

print("Making genomes...")
genomes = [makeGenome(i) for i in range(ploy)]

print("Writing genomes to %s/genomes.fa..." % tdir)
with open(tdir + "/genomes.fa", "w") as outhandle:
    outhandle.write(">genomes_p%d_s%d_t%d\n" % (ploy, gs*chrs, ploy*gs*chrs))
    for i in genomes:
        outhandle.write("".join(i))
print("Genome done.\n")

