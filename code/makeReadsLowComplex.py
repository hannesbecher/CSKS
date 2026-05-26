# make sequencing reads from genomes file (generated using makeGenomes.py)

# USAGE: python makeReadsLowComplex.py tdir cplx [--read-length 150] [--depth 50] [--error-prob 0.0 | --error-profile profile.txt]
import argparse

parser = argparse.ArgumentParser(
    description="Make low-complexity sequencing reads from a genomes.fa file."
)
parser.add_argument("tdir")
parser.add_argument(
    "cplx",
    type=float,
    help="Complexity of seq library, proportion of bases that have reads starting.",
)
parser.add_argument(
    "--read-length",
    "--rlen",
    dest="rlen",
    type=int,
    default=150,
    help="Read length (default: 150).",
)
parser.add_argument(
    "--depth",
    "--coverage",
    dest="cov",
    type=int,
    default=50,
    help="Sequencing coverage per haploid genome (default: 50).",
)
error_group = parser.add_mutually_exclusive_group()
error_group.add_argument(
    "--error-prob",
    dest="error_prob",
    type=float,
    default=0.0,
    help="Fixed per-nucleotide sequencing error probability (default: 0.0).",
)
error_group.add_argument(
    "--error-profile",
    dest="error_profile",
    help="Text file with one per-nucleotide error probability for each read position.",
)
args = parser.parse_args()

tdir = args.tdir
cplx = args.cplx
cov = args.cov
rlen = args.rlen
error_prob = args.error_prob
error_profile = args.error_profile

if not 0 <= cplx <= 1:
    parser.error("cplx must be between 0.0 and 1.0")
if cov <= 0:
    parser.error("depth must be a positive integer")
if rlen <= 0:
    parser.error("read length must be a positive integer")
if not 0 <= error_prob <= 1:
    parser.error("error probability must be between 0.0 and 1.0")

import numpy as np


def format_probability(value):
    return "{:.6g}".format(value)


def load_error_profile(path, expected_length):
    with open(path, "r") as infile:
        fields = infile.read().replace(",", " ").split()

    try:
        profile = [float(value) for value in fields]
    except ValueError:
        parser.error("error profile must contain only numeric probabilities")

    if len(profile) != expected_length:
        parser.error(
            "error profile length (%d) must match read length (%d)"
            % (len(profile), expected_length)
        )
    if any(prob < 0 or prob > 1 for prob in profile):
        parser.error("error profile probabilities must be between 0.0 and 1.0")

    return np.array(profile)


BASES = np.array(list("ACGT"))
ALT_BASES = {
    "A": np.array(list("CGT")),
    "C": np.array(list("AGT")),
    "G": np.array(list("ACT")),
    "T": np.array(list("ACG")),
}


def add_sequencing_errors(seq, error_probs):
    error_sites = np.random.random(len(seq)) < error_probs
    if not error_sites.any():
        return seq

    seq = list(seq)
    for pos in np.flatnonzero(error_sites):
        base = seq[pos].upper()
        seq[pos] = np.random.choice(ALT_BASES.get(base, BASES))
    return "".join(seq)


if error_profile:
    error_probs = load_error_profile(error_profile, rlen)
    error_mode = "profile"
else:
    error_probs = np.full(rlen, error_prob)
    error_mode = "fixed"
add_errors = error_probs.any()

print("########################################\nREAD MAKING SCRIPT")
print("Input parameters:")
print("  Temporary directory: %s" % tdir)
print("  Library complexity: %s" % cplx)
print("  Read length: %d" % rlen)
print("  Sequencing depth: %d" % cov)
print("  Error mode: %s" % error_mode)
if error_profile:
    print("  Error profile: %s" % error_profile)
    print("  Error profile positions: %d" % len(error_probs))
    print("  Error profile min probability: %s" % format_probability(error_probs.min()))
    print("  Error profile max probability: %s" % format_probability(error_probs.max()))
    print("  Error profile mean probability: %s" % format_probability(error_probs.mean()))
else:
    print("  Error probability: %s" % format_probability(error_prob))

print("Reading genomes...")
# the file has two lines: header and data
with open(tdir + "/genomes.fa", "r") as infile:
    hd, dat = infile.readlines() # there are only two lines: header and data
dat = dat.rstrip("\n")

pl = int(hd[:-1].split("_")[1][1:])    # get ploidy from header line
hgs = int(hd[:-1].split("_")[2][1:])   # get haploid genome size from header line
tlen = int(hd[:-1].split("_")[-1][1:]) # get total genome size from header line

# compute number of reads to generate
nreads = hgs * cov * pl // rlen
print("Genome parameters:")
print("  Ploidy: %d" % pl)
print("  Haploid genome size: %d" % hgs)
print("  Total genome size: %d" % tlen)
print("  Reads to generate: %d" % nreads)

# sample fragemtn within range(tlen-rlen) with probability cplx, without replacement
if cplx < 1.0:
    libStarts = np.random.choice(range(tlen-rlen), int(cplx*tlen), replace=False)
else:
    libStarts = range(tlen-rlen)
print("  Candidate read starts: %d" % len(libStarts))

# sample read locations
rstarts = np.random.choice(libStarts, nreads, replace=True)

print("Writing reads to %s/reads.fa..." % tdir)
with open(tdir + "/reads.fa", "w") as outfile:
    c = 0
    for i in rstarts:
        read = dat[i:(i+rlen)]
        if add_errors:
            read = add_sequencing_errors(read, error_probs)
        #outfile.write(">%d\n%s\n" % (c, dat[i:(i+rlen+1)]))
        outfile.write(">%d\n%s\n" % (c, read))
        c += 1
  
print("Reads done.\n")
