# count all paired mers in a gapped window

# can even have a weighted matrix for more important proteins in heptamer
# for example the CAC in the heptamer can be weighted

using BioSequences

export initialize_kmers
export as_index
export kmer_pair_count
export kmer_pair_count!

function initialize_kmers(seq::DnaSeq, k::Int, Nt_bits::DnaBits = NUCLEOTIDE_BITS)
    kmer = unsigned(0)
    for nt in view(seq, 1:k-1)
        kmer = (kmer << 2) + Nt_bits[nt]
    end
    return kmer, kmer
end

function as_index(kmer1::UInt, kmer2::UInt, k::Int)
    return ((kmer2 << (k << 1)) | kmer1) + 1
end

"""
    kmer_pair_count(seq::DnaSeq, k::Int = 3, Nt_bits::DnaBits = NUCLEOTIDE_BITS)

Counts all paired kmers in a Biosequences DNA longsequence (or a longsubseq),
where each index of the resulting `paired kmer frequency vector` correspond to the first and second kmers appended together.

For example, the 2mer pair `AT` and `GC` in a sequence `ATGC` would have a count of 1 in the resulting vector of length 4^(2*2) = 1
It would be at the index in the vector equivalent to `as_UInt(dna"ATGC")`
"""
function kmer_pair_count(seq::DnaSeq, k::Int = 3, Nt_bits::DnaBits = NUCLEOTIDE_BITS)
    n, bins, mask = length(seq), zeros(4^(2k)), unsigned((4^k)-1)
    kmer_i::UInt, kmer_j::UInt = initialize_kmers(seq, k)
    seqview = view(seq, k:n)

    # iterate to count pair mers
    for nt_i in seqview
        kmer_i = ((kmer_i << 2) & mask) + Nt_bits[nt_i]
        for nt_j in seqview
            kmer_j = ((kmer_j << 2) & mask) + Nt_bits[nt_j]
            bins[as_index(kmer_i, kmer_j, k)] += 1
        end
    end
    return bins
end

function kmer_pair_count!(seq::DnaSeq, k::Int, n::Int,
    bins::Vector{Float64}, mask::UInt, Nt_bits::DnaBits = NUCLEOTIDE_BITS)

    kmer_i::UInt, kmer_j::UInt = initialize_kmers(seq, k)
    seqview = view(seq, k:n)

    for nt_i in seqview
        kmer_i = ((kmer_i << 2) & mask) + Nt_bits[nt_i]
        for nt_j in seqview
            kmer_j = ((kmer_j << 2) & mask) + Nt_bits[nt_j]
            bins[as_index(kmer_i, kmer_j, k)] += 1
        end
    end
end

# interestingly, due to the sparsity it might (but probably not) be wise to not use bins