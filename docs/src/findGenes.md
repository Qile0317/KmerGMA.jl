# Homology searching
The main function of the package. Based on a novel kmer-based algorithm to reduce runtime over a whole genome to be *O(#bp)*.

The algorithm intakes a set of similar reference sequences - for example a collection of V-genes - and iterates over a genome assembly(ies) to find genes that are homologous to the references. 

```@docs
KmerGMA.findGenes
```

As mentioned in the docstring, if there surely are no "N" nucleotides present, setting `HasN` to `false` runs an alternate, much faster algorithm.

The `thr` argument can be toyed around with, as increasing the kmer distance threshold doesn't nessecarily always increase the amount of matches. The default distance is approximate against random sequences and may not always be the best metric.

!!! note
    The unoptimized algorithm that accounts for N nucleotides runs in just under 5 minutes to iterate over a genome of almost 4 billion bps. More benchmarking is on the way. More information is upcoming in a pre-print.