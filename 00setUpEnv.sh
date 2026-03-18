#conda create -n CSKS -c bioconda -c conda-forge msprime kmc python=3 numpy matplotlib wgsim ipython bwa-mem2 picard mosdepth samtools
#conda create -n CSKS_nomap -c bioconda -c conda-forge msprime kmc python=3 numpy matplotlib wgsim ipython
micromamba create -n CSKS -c conda-forge -c bioconda msprime=1.4 kmc python=3 numpy matplotlib wgsim ipython bwa-mem2 picard mosdepth samtools