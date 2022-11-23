# KmerGMA
[![Build Status](https://github.com/Qile0317/KmerGMA.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/Qile0317/KmerGMA.jl/actions/workflows/CI.yml?query=branch%3Amaster) [![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/Qile0317/LICENSE) [![Latest Release](https://img.shields.io/github/release/Qile0317/KmerGMA.jl.svg)](https://github.com/Qile0317/KmerGMA.jl/releases/latest) [![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://qile0317.github.io/KmerGMA.jl/dev/)

*A k-mer based approach for locating genes from gene families*

A package for finding homologues of a reference gene family/sequence using kmer manipulation swiftly. Also includes utilities such as an algorithmic to find all exact matches of a reference sequence to a genome. 

note:
- The package is in a very preliminary stage and can be heavily optimized. For potential developers, a pre-print of the main kmer-based algorithm will be available very soon. 

## Package Features
- Find approximate matches from a gene fammily (such as V(D)J genes) in a genome of a set of (or 1) reference sequence in linear time. 
- Find exact gene matches in a genome from a query sequence
- Utilities for kmer-tricks

## Installation
The package is not yet registered in the pkg registry. Until then, please use:

```julia
using Pkg
Pkg.add(PackageSpec(name="KmerGMA", url = "https://https://github.com/Qile0317/KmerGMA.jl.git"))
```

If you are interested in the cutting edge of the development, please check out the master branch to try new features before release.

## Documentation
Documentation is deployed at https://qile0317.github.io/KmerGMA.jl/dev/

## Performance
The current version of the approximate gene matching algorithm in the package processed a whole alpaca genome in just under 8 minutes over 4 billion base pairs. An alternate, unoptimized version ran in just under 5 minutes and a proposed, optimal algorithm has the potential to run in under a minute. More details will be shown very soon in a pre-print. 

## Contributions
The work was initially begun at Karolinska Institutet, under an in-progress project ```Discovery of Novel Germline Immunoglobulin alleles``` in which 2 approaches were utilized to expand camelid V(D)J databases. Thanks to @murrellb for massive support. More information is found at https://github.com/Qile0317/SoFoCompBio22