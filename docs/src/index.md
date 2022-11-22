# KmerGMA.jl
[![Build Status](https://github.com/Qile0317/KmerGMA.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/Qile0317/KmerGMA.jl/actions/workflows/CI.yml?query=branch%3Amaster) [![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/Qile0317/LICENSE)

*A k-mer based approach for locating genes from gene families*

A work-in-progress package for finding homologues of a reference gene family/sequence using kmer manipulation swiftly. Also includes utilities such as an algorithmic to find all exact matches of a reference sequence to a genome. 

!!! note
    The package is in a very preliminary stage and can be heavily optimized. For potential developers, a pre-print of the main kmer-based algorithm will be available very soon. 

## Package Features
- Find approximate matches from a gene fammily (such as V(D)J genes) in a genome of a set of (or 1) reference sequence in linear time. 
- Find exact matches in a genome from a query sequence(s)
- Utilities for kmer-tricks

## Installation
The package is not yet registered in the pkg registry. Until then, please use:

```julia
using Pkg
Pkg.add(PackageSpec(name="KmerGMA", url = "https://https://github.com/Qile0317/KmerGMA.jl.git"))
```

If you are interested in the cutting edge of the development, please check out the master branch to try new features before release.