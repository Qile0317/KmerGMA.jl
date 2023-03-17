using BioSequences
import Core

"""
    kmer_count(str::DnaSeq, k::Int, Nt_bits::DnaBits = NUCLEOTIDE_BITS)

Simple kmer counting function for a DNA BioSequence `str`, where `k` is the kmer length to count.
The function returns a kmer frequency vector where each INDEX of the vector corresponds to a unique kmer.
For each index, the unsigned bits of the index correspond to a binary representation of sequences where two bits encode one nucleotide.
To see what kmer each index corresponds to, see the function `as_kmer`
"""
function kmer_count(str::DnaSeq, k::Int, Nt_bits::DnaBits = NUCLEOTIDE_BITS)
    bins, mask, kmer = zeros(4^k), unsigned((4^k)-1), unsigned(0)
    for c in view(str, 1:k-1)
        kmer = (kmer << 2) + Nt_bits[c]
    end
    for c in view(str, k:length(str))
        kmer = ((kmer << 2) & mask) + Nt_bits[c]
        bins[kmer + 1] += 1
    end
    return bins
end

export kmer_count

# kmer counter that mutates the parameters
function kmer_count!(; str::DnaSeq, k::Int, 
    bins::Vector, mask::UInt,
    Nt_bits::DnaBits = NUCLEOTIDE_BITS)

    kmer = unsigned(0)
    for c in view(str, 1:k-1)
        kmer = (kmer << 2) + Nt_bits[c]
    end
    for c in view(str, k:length(str))
        kmer = ((kmer << 2) & mask) + Nt_bits[c]
        bins[kmer + 1] += 1
    end
end

export kmer_count!

# pairwise kmer distance of two sequences, not optimized for performance
function kmer_dist(seq1, seq2, k::Int, Nt_bits::DnaBits = NUCLEOTIDE_BITS)
    return (1/(2*k))*Distances.sqeuclidean(kmer_count(seq1, k, Nt_bits), kmer_count(seq2, k, Nt_bits))
end

function kmer_dist(seq1::DnaSeq, KFV::Kfv, k::Int, Nt_bits::DnaBits = NUCLEOTIDE_BITS)
    return (1/(2*k))*Distances.sqeuclidean(kmer_count(seq1, k, Nt_bits), KFV)
end

export kmer_dist

# code below are for user convenience with kmer and bit conversions

const BitNtDict = Dict{UInt, Seq}(
    unsigned(0) => dna"A", 
    unsigned(2) => dna"C", # intentional index switch 
    unsigned(1) => dna"G",
    unsigned(3) => dna"T")

"""
    get_kmer(kmer_uint::UInt, kmer_len::Int,Nt_bits::Dict{UInt, Seq} = BitNtDict)

Takes a unsigned integer representing a 2-bit-based dna kmer, and the actual length of the kmer, and returns the BioSequences dna sequence of the kmer. 
The Nt_bits argument can be ignored. Additionally, note that both input parameters will be modified to 0 or 1.
"""
function as_kmer(kmer_uint::UInt, kmer_len::Int, Nt_bits::Dict{UInt, Seq} = BitNtDict)
    output_seq = dna""
    while kmer_len > 0 ; curr_bit_nt = unsigned(0)
        for _ in 1:2
            remainder = kmer_uint % 2
            curr_bit_nt = (curr_bit_nt << 1) | remainder
            kmer_uint = unsigned(Int((kmer_uint-remainder)*0.5))
        end
        output_seq = Nt_bits[curr_bit_nt] * output_seq
        kmer_len -= 1
    end
    return output_seq
end

export as_kmer

"""
    UInt(kmer_seq::LongSequence{DNAAlphabet{4}}, Nt_bits::DnaBits = NUCLEOTIDE_BITS)

Convert a dna biosequence `kmer_seq` into an unsigned integer. Users can ignore the second argument
"""
function as_UInt(kmer_seq::Seq, Nt_bits::DnaBits = NUCLEOTIDE_BITS)
    bit_kmer = unsigned(0)
    for c in kmer_seq
        bit_kmer = ((bit_kmer << 2) | Nt_bits[c])
    end
    return bit_kmer
end

export as_UInt