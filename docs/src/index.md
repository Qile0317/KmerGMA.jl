# KmerGMA.jl
*A k-mer based approach for locating genes from gene families*

## Package Features
- Find exact matches in a genome from a query sequence
- Find approximate matches in a genome of a set of (or 1) reference sequence.

## Function Documentation
```@docs
    findgenes()
```
The only relevant function for users is ```findgenes()``` which takes a reference set of/individual sequence(s) and scans a genome for approximate homologous matches. This can be applied for finding V genes, for example. 
