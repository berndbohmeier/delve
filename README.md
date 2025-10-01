[![CI](https://github.com/berndbohmeier/delve/actions/workflows/ci.yml/badge.svg)](https://github.com/berndbohmeier/delve/actions/workflows/ci.yml)
[![Conda Downloads](https://img.shields.io/conda/dn/bioconda/delve-bio?color=green&link=https%3A%2F%2Fanaconda.org%2Fbioconda%2Fdelve-bio)](https://anaconda.org/bioconda/delve-bio)
![OS - Linux | OSX](https://img.shields.io/badge/OS-Linux_|_OSX-informational)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/delve-bio/README.html)

# delve
Variant caller for mixed infections

## Install
Delve can be install from bioconda

```
conda install bioconda::delve-bio
```

## Run
The minimum you need to run delve is an indexed BAM file with aligned reads to a reference genome, for example with minimap2, a reference genome in fasta format with index. Optionally you can provide a BED file with regions over which to call the variants.

```
delve -f reference.fasta [-R amplicons.bed] input.bam
```

## Features
Delve internally uses [rust-htslib](https://github.com/rust-bio/rust-htslib) to build a pileup of reads over each position, then uses a statistical model, which uses the base read errors of the reference and alternative bases on this position to calcululate a likelihood of a proportion p of alt vs. ref bases, and then uses a [Likelihood-ratio Test](https://en.wikipedia.org/wiki/Likelihood-ratio_test) to determine if there is enough evidence that a variant exists.
After, it uses a strand bias test to filtered out wrong calls.
Currently delve only calls variants on one sample at a time.
