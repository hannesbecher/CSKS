# Python script to make a diploid genome with a certain heterozygosity

# USAGE: python makeDiploidNoCoal.py gs theta gc tdir propDel
# TRY: python makeDiploidNoCoal.py 100000 0.02 0.5 abc123 0.1
import numpy as np
import random
import sys
import copy

gs = int(sys.argv[1])   # size of one chromosome
theta = float(sys.argv[2])
gc = float(sys.argv[3])
tdir = sys.argv[4]
pDel = float(sys.argv[5])
chrs = 1

print("Making random genome...")
rng = np.random.default_rng()
refArr = np.random.choice(["A","T","C","G"], gs, replace=True, p=[(1-gc)/2,(1-gc)/2, (gc)/2, (gc)/2 ])


def pickAlt(ref, gc=gc):
    a=rng.choice(["A","T","C","G"],1,replace=False,p=[(1-gc)/2,(1-gc)/2, (gc)/2, (gc)/2 ])[0]
    if a == ref:
        a = pickAlt(ref, gc=gc)
    return a

print("Make heterozygous sites...")    
numHet = np.random.binomial(gs, theta)
#hetPos = random.sample(range(gs), numHet)
hetPos = list(set(np.random.choice(range(gs), numHet, replace=True)))
altArr = copy.deepcopy(refArr)
for i in hetPos:
    altArr[i] = pickAlt(refArr[i])


print("Generating deletions...")    
# just one for a start
dNt = [int(gs * pDel)] * chrs # number of nts in the deletion (for each chromosome)
print(dNt)
dStarts = [random.sample(range(gs-dNt[i]), 1)[0] for i in range(chrs)]
altArrList = [altArr]
#altArrList = [np.hstack([
#     altArrList[i][0:(dStarts[i])],
#     altArrList[i][(dStarts[i]+dNt[i]):-1]
#     ]) for i in range(chrs)]
# for 1 chr only
altArrList = [np.delete(altArrList[0], np.arange(dStarts[0], dStarts[0]+ dNt[0]))]




print("Writing genomes...")
with open(tdir + "/genomes.fa", "w") as outhandle:
    outhandle.write(">genomes_p%d_s%d_t%d\n" % (2, gs*chrs, 2*gs*chrs))
    outhandle.write("".join(refArr))
    #outhandle.write("".join(altArr))
    outhandle.write("".join(altArrList[0]))
with open(tdir + "/reference.fa", "w") as outhandle:
    outhandle.write(">ref\n")
    outhandle.write("".join(refArr))
print("Genome done.\n")
    
