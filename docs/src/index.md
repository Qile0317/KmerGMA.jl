# KmerGMA
[![Build Status](https://github.com/Qile0317/KmerGMA.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/Qile0317/KmerGMA.jl/actions/workflows/CI.yml?query=branch%3Amaster) [![Codecov test coverage](https://codecov.io/gh/Qile0317/KmerGMA.jl/branch/master/graph/badge.svg)](https://app.codecov.io/gh/Qile0317/KmerGMA.jl?branch=master) [![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/Qile0317/KmerGMA.jl/LICENSE) [![Latest Release](https://img.shields.io/github/release/Qile0317/KmerGMA.jl.svg)](https://github.com/Qile0317/KmerGMA.jl/releases/latest) [![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://qile0317.github.io/KmerGMA.jl)

*A k-mer based approach for homolgy searching*

A package for finding homologues of a reference gene family/sequence using kmer manipulation with a single scan. Also includes utilities such as an algorithm to find all exact matches of a reference sequence to a genome. 

The approach was developed to predict novel VDJ alelles in genomes based on sets of reference sequences, but should be just as viable for other gene homologues that arise from gene families.

## Package Features
- Homology searching of 1 query sequence or a collection of queries to a local database(s)
- Find exact gene matches in a genome from a query sequence
- Many subsequence counting utilities for kmers, paired kmers, and strobemers (https://doi.org/10.1101/gr.275648.121)

## Installation
The package is not yet registered in the pkg registry. Until then, please use:

```julia
using Pkg
Pkg.add(PackageSpec(name="KmerGMA", url = "https://github.com/Qile0317/KmerGMA.jl.git"))
```

If you are interested in the cutting edge of the development, please check out the master branch to try new features before release.

## Usage Example
To conduct homology searching of a set of sequences in a local fasta file with another query sequence fasta file, simply do:
```julia
using KmerGMA
KmerGMA.findGenes(genome_path = "my_sequences.fasta", ref_path = "my_query_sequence_family.fasta")
```
Where `genome_path` is the file location of a fasta file containing sequences to search over, and `ref_path` is the file location of a fasta file containing the query sequence or a set of alike query sequences (For example V-genes). 

The function defaults to returning a vector of fasta records. 

Alternatively, if accuracy is favored over speed, then its suggested to be running the following: 
```julia
KmerGMA.findGenes_cluster_mode(genome_path = "my_sequences.fasta", ref_path = "my_query_sequence_family.fasta")
```
Where speed is sacrificed for accuracy.

The functions have many optional paramaters to optimize/adjust its performance. See the documentation for more details.

## Performance
The current version of the homology searching `findGenes` function can iterate on average 40 megabases per second. So it would take about 80 seconds for the human genome. The performance of `findGenes_cluster_mode` slows porportionally to the number of reference sequence clusters, so for 5 clusters it would be 40/5 = 8 megabases/second. 

## Contributions
The work was initially begun at Karolinska Institutet, as a side project to the in-progress project ```Discovery of Novel Germline Immunoglobulin alleles``` in which 2 approaches were utilized to expand camelid V(D)J databases. Thanks to @murrellb for massive support. More information is found at https://github.com/Qile0317/SoFoCompBio22