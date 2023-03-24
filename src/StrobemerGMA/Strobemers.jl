# implementation of RandStrobes
# This is a super slow and unoptimized MWE, future versions will be bit based
# I think the big thing here is that strobemers give a simple way to map one kmer to another that skips indels hopefully

using Plots, Random, Distances

include("Consts.jl")
include("Kmers.jl")

function randstrobe_score(s1::DnaSeq, s2::DnaSeq, q::Int)
    return (as_UInt(s1) + as_UInt(s2)) % q
end

# can be generalized
function get_strobe_2_mer(
    seq::DnaSeq, s::Int = 2,
    w_min::Int = 3, w_max::Int = 5,
    q::Int = 5, # q should be prime
    ; withGap::Bool = false) 

    first_strobe = Seq(seq[1:s])
    min_score::Int = 2 << 63
    min_ind::Int = w_min
    for i in w_min:w_max
        curr_score = randstrobe_score(first_strobe, view(seq, i:i+s-1), q)
        if curr_score <= min_score # do <= so that strobes get further! thats usually good
            min_score = curr_score
            min_ind = i 
        end
    end
    if !withGap
        return first_strobe * Seq(view(seq, min_ind:min_ind+s-1))
    end
    return first_strobe * dna"-"^(min_ind-s-1) * Seq(seq[min_ind:min_ind+s-1]) * dna"-"^(length(seq)-min_ind-s+1)
end

# initial strobemer count
function strobe_2_mer_count(
    seq::DnaSeq, s::Int = 2,
    w_min::Int = 3, w_max::Int = 5,
    q::Int = 5)

    k = w_max+s-1
    bins = zeros(4^k)
    for i in 1:length(seq)-k+1
        strobemer = get_strobe_2_mer(view(seq,i:i+k-1), s, w_min, w_max, q)
        bins[as_UInt(strobemer) + 1] += 1
    end
    return bins
end

function dist(kfv1::Vector{Float64}, kfv2::Vector{Float64})
    return (1/(2*(log(4, length(kfv1)))))*Distances.sqeuclidean(kfv1, kfv2)
end

"""
Random.seed!(42)
vgene = dna"caggtccagctggtgcagccaggggctgagctgaggaagcctggggctttgctgaaggtctcctgcaaggcttctggatacaccttcaccagctactacatagactgggtgcgacaggcccctggacaagggcttgggtgggtgggaagaattgaccctgaagacggtggcacaaattatgcacagaagttccagggcagagtcaccttgactgcagacacgtccaccagcacagcctacgtggagctgagcagtctgagatctgaggacacggccgtgtgttactgtgtgagaga"
vgene_10_div = mutate_seq(vgene, 0.3)
rand_seq = randdnaseq(296)

seq1 = strobe_2_mer_count(vgene)
seq2 = strobe_2_mer_count(vgene_10_div)
rand_seq3 = strobe_2_mer_count(rand_seq)
plot(seq1); plot!(seq2)
dist(seq1,seq2)
dist(seq1, rand_seq3)
"""