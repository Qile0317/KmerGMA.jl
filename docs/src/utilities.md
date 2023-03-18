# Utility functions for processing sequences and kmers
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
```

## paired kmer utilities
The following function count all kmer pairs within a sequence, which is a feature to be further utilized in an upcoming release. Counting paired kmers could serve as a more robust alternative to regular kmers in many applications.
```@docs
KmerGMA.kmer_pair_count
```