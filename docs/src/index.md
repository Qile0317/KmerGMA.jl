# KmerGMA.jl
[![Build Status](https://github.com/Qile0317/KmerGMA.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/Qile0317/KmerGMA.jl/actions/workflows/CI.yml?query=branch%3Amaster) [![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/Qile0317/LICENSE) [![Latest Release](https://img.shields.io/github/release/Qile0317/KmerGMA.jl.svg)](https://github.com/Qile0317/KmerGMA.jl/releases/latest) [![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://qile0317.github.io/KmerGMA.jl/)

*A k-mer based approach for locating genes from gene families*

A work-in-progress package for finding homologues of a reference gene family/sequence using kmer manipulation swiftly. Also includes utilities such as an algorithmic to find all exact matches of a reference sequence to a genome. 

!!! note
    The package is in a very preliminary stage and can be heavily optimized. For potential developers, a pre-print of the main kmer-based algorithm will be available very soon. 

## Package Features
- Homology searching of 1 query sequence or a collection of queries to a local database(s)
- Find exact gene matches in a genome from a query sequence
- Utilities for kmer-tricks

## Installation
The package is not yet registered in the pkg registry. Until then, please use:

```julia
using Pkg
Pkg.add(PackageSpec(name="KmerGMA", url = "https://github.com/Qile0317/KmerGMA.jl.git"))
```

If you are interested in the cutting edge of the development, please check out the master branch to try new features before release.

## Performance
The current version of the approximate gene matching algorithm in the package can iterate over a collection of sequences of around 4 billion basepairs in 5 minutes. Current benchmarking indicates that optimziations can improve this by a full minute in subsequent updates. Proposed optimizations has the theoretical potential to run in under a minute. More details will be shown very soon in a pre-print. 