# Homology searching
The main function of the package. Based on a novel kmer-based algorithm to reduce runtime over a whole genome to be $O(#bp)$.

The algorithm intakes a set of similar reference sequences - for example a collection of V-genes - and iterates over a genome assembly(ies) to find genes that are homologous to the references. More information is upcoming in a pre-print.

```@docs
KmerGMA.findGenes
```
!!! note
    The unoptimized algorithm runs in just under 8 minutes to iterate over a genome of almost 4 billion bps. More benchmarking is on the way.