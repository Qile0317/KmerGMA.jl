using BioSequences

function kmer_count(str::DnaSeq, k::Int, Nt_bits::DnaBits = NUCLEOTIDE_BITS)
    bins, mask, kmer = zeros(4^k), unsigned((4^k)-1), unsigned(0)
    for c in str[1:k-1]
        kmer = (kmer << 2) + Nt_bits[c]
    end
    for c in str[k:end]
        kmer = ((kmer << 2) & mask) + Nt_bits[c]
        bins[kmer + 1] += 1
    end
    return bins
end

function kmer_count!(; str::DnaSeq, k::Int, 
    bins::Vector, mask::UInt,
    Nt_bits::DnaBits = NUCLEOTIDE_BITS)

    kmer = unsigned(0)
    for c in str[1:k-1]
        kmer = (kmer << 2) + Nt_bits[c]
    end
    for c in str[k:end]
        kmer = ((kmer << 2) & mask) + Nt_bits[c]
        bins[kmer + 1] += 1
    end
end

# pairwise kmer distance of two sequences, not optimized for performance
function kmer_dist(seq1, seq2, k::Int, Nt_bits::DnaBits = NUCLEOTIDE_BITS)
    return (1/(2*k))*Distances.sqeuclidean(kmer_count(seq1, k, Nt_bits), kmer_count(seq2, k, Nt_bits))
end

function kmer_dist(seq1::DnaSeq, KFV::Kfv, k::Int, Nt_bits::DnaBits = NUCLEOTIDE_BITS)
    return (1/(2*k))*Distances.sqeuclidean(kmer_count(seq1, k, Nt_bits), KFV)
end