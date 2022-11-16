# KmerGMA.jl
*A k-mer based approach for locating genes from gene families*

## Package Features
- Find exact matches in a genome from a query sequence
- Find approximate matches in a genome of a set of (or 1) reference sequence.

## Function Documentation
```@docs
genKmers(k::Int64,
         Dictionary::Bool = true;
         withN::Bool = false,
         cumulative::Bool = false)
```
