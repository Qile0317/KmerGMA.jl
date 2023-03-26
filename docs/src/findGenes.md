# Homology searching
Here is the main function of the package. Based on a novel kmer-based algorithm to reduce runtime over a whole genome to be *O(#bp)* on average.

!!! note
    The documentation is unfinished and may be somewhat unclear. For more details on how the algorithm functions, some pre-print of the homology searching approach will be uploaded soon.

The algorithm intakes a set of similar reference sequences - for example a collection of V-genes - and iterates over a genome assembly(ies) to find genes that are homologous to the references. 

```@docs
KmerGMA.findGenes
```

An alternative version of `KmerGMA.findGenes` which clusters the reference sequence set into similar subsets is also avaliable, though its running speed would be "the number of subsets" times slower. However, the following function is recommended for when speed is less important than accuracy.

```@docs
KmerGMA.findGenes_cluster_mode
```

To write resulting vector of hits into a fasta file, use the following function
```@docs
KmerGMA.write_results
```

There is also possible pre-processing of the set of reference sequences that can be used in `kmerGMA.findGenes_cluster_mode` to further improve accuracy (though usually at the cost of a bit more more speed) More on this will be expanded on.