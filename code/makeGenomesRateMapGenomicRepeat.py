# run msprime to generate genome sequences
# alternative scipt, gnerate alternative genotypes for each site, not matter is variant or not.
# USAGE: python makeGenomes.py gs ploy theta
import msprime
import sys
import numpy as np
import random
import copy

print("########################################\nGENOME MAKING SCRIPT (WITH REC RATE MAP)")
rng = np.random.default_rng()

gc=0.5
gs = int(sys.argv[1])   # size of one chromosome (10 are simulated)!
ploy = int(sys.argv[2])
theta = float(sys.argv[3])
tdir = sys.argv[4]
recBaseRate = 10e-4

#tdir = "/tmp/tmp.0gVEKSPNrd"
#gs, ploy, theta = 1000000, 2, 0.005
#gs, ploy, theta = 1000000, 2, 0.01
print("Simulating data for haploid GS %d, ploidy %d, and theta %1.3f..." % (gs, ploy, theta))


chrs=10 # number of replicate runs corresponding to "chromosomes" of equal size

mapBreaks= [0]
for i in range(chrs):
    mapBreaks.append((i+1) * gs - 1)
    mapBreaks.append((i+1) * gs)
    
recRates = []
for i in range(chrs):
    recRates.append(recBaseRate)
    recRates.append(0.5)    
recMap = msprime.RateMap(position=mapBreaks, rate=recRates)

ts = msprime.sim_ancestry(
        samples=ploy,
        ploidy=1,
        population_size=1,
        discrete_genome=True,# default setting now?
        recombination_rate=recMap)#, # possibly adjust
#        sequence_length=gs) 
#        random_seed=123456)

print("Adding mutations...")


mutated_ts = [msprime.sim_mutations(i, rate=theta/2, model=msprime.JC69()) for i in [ts]]


# A list of lists, one for each "chromosome"
# Each chromsome list contains tuples of variant positions and genotypes
print("Extracting variant site information and per-site GTs...")
tupListList = [[(i.site.position, i.genotypes) for i in j.variants()] for j in mutated_ts]


varSites = dict()
for i in range(1):
    for j in tupListList[i]:
        varSites[int(j[0]) + (i * gs)] = np.random.choice(["A","C","G","T"], 4, replace=False)[j[1]]
print(len(varSites.keys()))

# 1D array of nucleotides:
print("Making a random genome ref, length %d..." % (gs*chrs))
refArr = np.random.choice(["A","C","G","T"], chrs*gs, replace=True, p=[(1-gc)/2,(1-gc)/2, (gc)/2, (gc)/2 ])

elementLen = 500
repElement = np.random.choice(["A","C","G","T"], elementLen, replace=True)

elemPos = np.random.randint(0, chrs*gs - elementLen, int(chrs * gs /2 / elementLen))
for i in elemPos:
    refArr[i:(i+elementLen)] = repElement

def makeGenome(k):
    a = copy.deepcopy(refArr)
    for i in varSites.keys():
        a[i] = varSites[i][k]
    return a

print("Making genomes...")
genomes = [makeGenome(i) for i in range(ploy)]

#sum([1 for i in range(10000000) if genomes[0][i] != genomes[1][i]])

print("Writing genomes to %s/genomes.fa..." % tdir)
with open(tdir + "/genomes.fa", "w") as outhandle:
    outhandle.write(">genomes_p%d_s%d_t%d\n" % (ploy, gs*chrs, ploy*gs*chrs))
    for i in genomes:
        outhandle.write("".join(i))
print("Genome done.\n")

