# run msprime to generate genome sequences
# USAGE: python makeGenomes.py gs ploy theta
import msprime
import sys
import numpy as np
import random
import copy
import tskit

print("########################################\nGENOME MAKING SCRIPT")
rng = np.random.default_rng()

gc=0.5
gs = int(sys.argv[1])   # size of one chromosome
ploy = int(sys.argv[2])
theta = float(sys.argv[3])
tdir = sys.argv[4]
pDel = float(sys.argv[5])

if ploy < 1:
    raise ValueError("ploidy must be at least 1")
if pDel < 0.0 or pDel > 1.0:
    raise ValueError("pDel must be between 0.0 and 1.0")
if pDel != 0.0 and ploy != 2:
    raise ValueError("Nonzero deletions are only supported for ploidy 2")

#tdir = "/tmp/tmp.0gVEKSPNrd"
#gs, ploy, theta = 10000, 2, 0.005
print("Simulating data for haploid GS %d, ploidy %d, and theta %1.3f..." % (gs, ploy, theta))

# one chromosome only for now
chrs=1 # number of replicate runs corresponding to "chromosomes" of equal size

ts = list(msprime.sim_ancestry(
        samples=ploy,
        ploidy=1,
        population_size=1,
        discrete_genome=True,# default setting now?
        recombination_rate=1e-3, # possibly adjust
        sequence_length=gs,
        num_replicates=chrs) )
#        random_seed=123456)

print("Adding mutations...")
# !! Here is a logical break all is made for only one chromosome from here on!!
mutated_ts = msprime.sim_mutations(ts[0], rate=theta/2)
#        random_seed=123456)


asFasta = mutated_ts.as_fasta(reference_sequence=tskit.random_nucleotides(mutated_ts.sequence_length),
                              wrap_width=0,
                              isolated_as_missing=False).split("\n")
haplotypes = asFasta[1::2]
if len(haplotypes) != ploy:
    raise RuntimeError("Expected %d haplotypes, found %d" % (ploy, len(haplotypes)))

# make one deletion, at end of chr
dLen = int(pDel * gs)
total_length = ploy * gs * chrs

print("Writing genomes to %s/genomes.fa..." % tdir)
with open(tdir + "/genomes.fa", "w") as outhandle:
    if dLen != 0:
        haplotypes[1] = haplotypes[1][:-dLen]
        total_length -= dLen
    genome_sequence = "".join(haplotypes)
    if len(genome_sequence) != total_length:
        raise RuntimeError("Genome sequence length does not match header length")
    outhandle.write(">genomes_p%d_s%d_t%d\n" % (ploy, gs*chrs, total_length))
    outhandle.write(genome_sequence + "\n")
print("Genome done.\n")
