# make sequencing reads from genomes file (generated using makeGenomes.py),
# inserting an "imperfect PCR" amplification step between library-complexity
# sampling and read sampling.
#
# Model: the library consists of `cplx * tlen` candidate template molecules
# (the same "candidate read starts" concept as makeReadsLowComplex.py). In
# the "perfect PCR" model (makeReadsLowComplex.py) every molecule is equally
# available for read sampling, so all duplication comes purely from the
# finite pool of molecules (cplx < 1). Here, before read sampling, each
# molecule is additionally given an amplification factor drawn from a
# distribution (geometric to start with; log-normal, following Rochette et
# al. 2023, as a second option). Read sampling is then a weighted draw
# proportional to these amplification factors, so molecules amplified more
# heavily contribute proportionally more reads -- adding a second,
# efficiency-driven source of duplication on top of the library-complexity
# effect already present in the original model.
#
# USAGE: python makeReadsImperfectPCR.py tdir cplx [--read-length 150] [--depth 50]
#        [--error-prob 0.0 | --error-profile profile.txt]
#        [--pcr-dist none|geometric|lognormal] [--pcr-param P]
import argparse

parser = argparse.ArgumentParser(
    description="Make low-complexity sequencing reads with an imperfect-PCR "
    "amplification step from a genomes.fa file."
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
parser.add_argument(
    "--pcr-dist",
    dest="pcr_dist",
    choices=["none", "geometric", "lognormal"],
    default="none",
    help="Distribution of per-molecule PCR amplification factors (default: none, "
    "i.e. the original 'perfect PCR' behaviour of makeReadsLowComplex.py).",
)
parser.add_argument(
    "--pcr-param",
    dest="pcr_param",
    type=float,
    default=None,
    help="Parameter of the amplification-factor distribution. For "
    "--pcr-dist geometric this is the per-trial success probability p in "
    "(0, 1] (mean copies = 1/p; p=1 means no amplification noise, i.e. "
    "equivalent to --pcr-dist none). For --pcr-dist lognormal this is the "
    "log-scale standard deviation sigma >= 0 (sigma=0 means no amplification "
    "noise). Required unless --pcr-dist none.",
)
args = parser.parse_args()

tdir = args.tdir
cplx = args.cplx
cov = args.cov
rlen = args.rlen
error_prob = args.error_prob
error_profile = args.error_profile
pcr_dist = args.pcr_dist
pcr_param = args.pcr_param

if not 0 <= cplx <= 1:
    parser.error("cplx must be between 0.0 and 1.0")
if cov <= 0:
    parser.error("depth must be a positive integer")
if rlen <= 0:
    parser.error("read length must be a positive integer")
if not 0 <= error_prob <= 1:
    parser.error("error probability must be between 0.0 and 1.0")

if pcr_dist == "none":
    if pcr_param is not None:
        parser.error("--pcr-param must not be set when --pcr-dist is none")
elif pcr_dist == "geometric":
    if pcr_param is None:
        parser.error("--pcr-param (success probability p, 0 < p <= 1) is required for --pcr-dist geometric")
    if not 0 < pcr_param <= 1:
        parser.error("for --pcr-dist geometric, --pcr-param (p) must be in (0, 1]")
elif pcr_dist == "lognormal":
    if pcr_param is None:
        parser.error("--pcr-param (sigma, sigma >= 0) is required for --pcr-dist lognormal")
    if pcr_param < 0:
        parser.error("for --pcr-dist lognormal, --pcr-param (sigma) must be >= 0")

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

print("########################################\nREAD MAKING SCRIPT (imperfect PCR)")
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
print("  PCR amplification distribution: %s" % pcr_dist)
if pcr_dist != "none":
    print("  PCR amplification parameter: %s" % format_probability(pcr_param))

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

# sample fragment start positions within range(tlen-rlen) with probability cplx,
# without replacement -- these are the candidate template molecules
if cplx < 1.0:
    libStarts = np.random.choice(np.arange(tlen - rlen), int(cplx * tlen), replace=False)
else:
    libStarts = np.arange(tlen - rlen)
n_molecules = len(libStarts)
print("  Candidate read starts (template molecules): %d" % n_molecules)

# assign a per-molecule PCR amplification factor and turn it into sampling
# weights. pcr_dist == "none" reproduces the original uniform (perfect PCR)
# behaviour exactly (weights=None -> np.random.choice samples uniformly).
if pcr_dist == "none":
    weights = None
    amp = None
elif pcr_dist == "geometric":
    # np.random.geometric support is {1, 2, 3, ...}; mean = 1/p.
    # p=1 -> every molecule gets exactly 1 "copy", i.e. no amplification noise.
    amp = np.random.geometric(pcr_param, size=n_molecules).astype(float)
    weights = amp / amp.sum()
elif pcr_dist == "lognormal":
    # sigma=0 -> every molecule gets amplification factor 1, i.e. no
    # amplification noise (mean of the underlying normal is fixed at 0 so the
    # median amplification factor is always 1 regardless of sigma).
    amp = np.random.lognormal(mean=0.0, sigma=pcr_param, size=n_molecules)
    weights = amp / amp.sum()

if amp is not None:
    print("  Amplification factor summary: mean=%.4g, sd=%.4g, min=%.4g, max=%.4g"
          % (amp.mean(), amp.std(), amp.min(), amp.max()))
    print("  Amplification factor CV: %.4g" % (amp.std() / amp.mean()))

# sample read locations (weighted by amplification factor when pcr_dist != none)
rstarts = np.random.choice(libStarts, nreads, replace=True, p=weights)

print("Writing reads to %s/reads.fa..." % tdir)
with open(tdir + "/reads.fa", "w") as outfile:
    c = 0
    for i in rstarts:
        read = dat[i:(i+rlen)]
        if add_errors:
            read = add_sequencing_errors(read, error_probs)
        outfile.write(">%d\n%s\n" % (c, read))
        c += 1

print("Reads done.\n")
