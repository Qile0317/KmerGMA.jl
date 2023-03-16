# Homology searching
The main function of the package. Based on a novel kmer-based algorithm to reduce runtime over a whole genome to be *O(#bp)*.

The algorithm intakes a set of similar reference sequences - for example a collection of V-genes - and iterates over a genome assembly(ies) to find genes that are homologous to the references. 

```@docs
KmerGMA.findGenes
KmerGMA.write_results
```

!!! note
    The documentation is unfinished and may be somewhat unclear. For more details on how the algorithm functions, some pre-print of the homology searching approach will be uploaded soon.