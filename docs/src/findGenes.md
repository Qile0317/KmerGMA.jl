# Homology searching
Here is the main function of the package. Based on a novel kmer-based algorithm to reduce runtime over a whole genome to be *O(#bp)* on average.

The algorithm intakes a set of similar reference sequences - for example a collection of V-genes - and iterates over a genome assembly(ies) to find genes that are homologous to the references. 

```@docs
KmerGMA.findGenes
```

!!! note
    The documentation is unfinished and may be somewhat unclear. For more details on how the algorithm functions, some pre-print of the homology searching approach will be uploaded soon.

To write the resulting hits into a fasta file, use the following function
```@docs
KmerGMA.write_results
```