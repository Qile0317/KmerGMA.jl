# Utility functions for processing sequences and subsequences
## fasta file utilities
```@docs
KmerGMA.fasta_id_to_cumulative_len_dict
KmerGMA.getSeq
```

## kmer utilities
```@docs
KmerGMA.kmer_count
KmerGMA.as_kmer
KmerGMA.as_UInt
KmerGMA.kmer_dist
```

## paired kmer utilities
The following function count all kmer pairs within a sequence, which is a feature to be further utilized in an upcoming release. Counting paired kmers could serve as a more robust alternative to regular kmers in many applications.
```@docs
KmerGMA.kmer_pair_count
```

## strobemer utilities (experimental)
Strobemers consist of two or more linked shorter k-mers, where the combination of linked k-mers is decided by a hash function. This implementation is currently being experimented with and is unoptimized, but may be implemented in the near future as a new homology searching algorithm. 

The strobemers implemented are simply randstrobes with two k-mers, and in the future will be a part of a new module `StrobemerGMA`
```@docs
KmerGMA.get_strobe_2_mer
KmerGMA.ungapped_strobe_2_mer_count
```

## sequence divergence (work in progress)
```@docs
KmerGMA.mutate_seq
KmerGMA.mutate_seq!
```