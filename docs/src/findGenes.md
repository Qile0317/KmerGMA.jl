# KmerGMA.findGenes
The star of the package. Based on a novel kmer-based algorithm to reduce runtime over a whole genome to be $O(#bp)$.

The algorithm intakes a set of similar reference sequences - for example a collection of V-genes - and looks at a genome assembly to find genes that are homologous to the references. More information is upcoming in a pre-print.

```@docs
KmerGMA.findGenes
```