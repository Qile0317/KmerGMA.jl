# implementation of RandStrobes
# This is a super slow and unoptimized MWE, future versions will be bit based
# I think the big thing here is that strobemers give a simple way to map one kmer to another that skips indels hopefully

using Plots, Random, Distances

export randstrobe_score
export get_strobe_2_mer
export ungapped_strobe_2_mer_count

# convinienice strobemer hashing function (should use uint and bitwise ops in the future)
function randstrobe_score(s1::DnaSeq, s2::DnaSeq, q::Int)
    return (as_UInt(s1) + as_UInt(s2)) % q
end

# an even simpler way is literally to just iterate and get the minimum number

# an alternative is to use syncmers (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7869670/)
# the minimizer can even be weighted (https://academic.oup.com/bioinformatics/article/36/Supplement_1/i111/5870473)
# maybe some analysis of common v gene stuff can be used

"""
    get_strobe_2_mer(
        seq::LongSequence{DNAAlphabet{4}},
        s::Int = 2,
        w_min::Int = 3,
        w_max::Int = 5,
        q::Int = 5;
        withGap::Bool = true) 

Get the randstrobe that consists of two kmers of the current sequence `seq`. 
The `withGap` named argument indicates whether to return the strobemer with
sequence gaps or just as a regular kmer with all gaps removed

...
# Strobemer parameters
- `s::Int = 2`: length of the kmers that consists of the strobemer.
- `w_min::Int = 3`: the start of the window to obtain the second kmer from.
- `w_max::Int = 5`: the end of the window to obtain the second kmer from.
- `q::Int = 5`: the prime number that the randstrobe hashing function should use. 
...

It is recommended to keep all optional strobemer parameters as is
"""
function get_strobe_2_mer(
    seq::DnaSeq, s::Int = 2,
    w_min::Int = 3, w_max::Int = 5,
    q::Int = 5;
    withGap::Bool = true) 

    first_strobe = Seq(seq[1:s])
    min_score::Int = 2 << 63
    min_ind::Int = w_min
    for i in w_min:w_max
        curr_score = randstrobe_score(first_strobe, view(seq, i:i+s-1), q)
        if curr_score <= min_score # do <= so that strobes get further!
            min_score = curr_score
            min_ind = i 
        end
    end
    if !withGap
        return first_strobe * Seq(view(seq, min_ind:min_ind+s-1))
    end
    return first_strobe * dna"-"^(min_ind-s-1) * Seq(seq[min_ind:min_ind+s-1]) * dna"-"^(length(seq)-min_ind-s+1)
end

"""
    ungapped_strobe_2_mer_count(
        seq::DnaSeq;
        s::Int = 2,
        w_min::Int = 3,
        w_max::Int = 5,
        q::Int = 5)

Get the ungapped randstrobe (two kmers of the current sequence `seq`) frequency vector. 
The vector is ``4^{2s}`` long, where each index corresponds to the strobemer with gaps removed.
So the strobemer `AC--GT--` would correspond to the index corresponding to `ACGT`. To convert
from index to the ungapped strobemer, see `KmerGMA.as_kmer`

...
# Strobemer parameters
- `s::Int = 2`: length of the kmers that consists of the strobemer.
- `w_min::Int = 3`: the start of the window to obtain the second kmer from.
- `w_max::Int = 5`: the end of the window to obtain the second kmer from.
- `q::Int = 5`: the prime number that the randstrobe hashing function should use. 
...

It is currently the responsibility of the user to ensure that all parameters do not conflict with eachother 
"""
function ungapped_strobe_2_mer_count(
    seq::DnaSeq; s::Int = 2,
    w_min::Int = 3, w_max::Int = 5,
    q::Int = 5)

    k = w_max+s-1
    bins = zeros(4^(2*s))
    for i in 1:length(seq)-k+1
        strobemer = get_strobe_2_mer(view(seq,i:i+k-1), s, w_min, w_max, q; withGap = false)
        bins[as_UInt(strobemer) + 1] += 1
    end
    return bins
end

function ungapped_strobe_2_mer_count!(
    seq::DnaSeq, bins::BinInput,
    s::Int = 2, w_min::Int = 3, w_max::Int = 5,
    q::Int = 5)

    k = w_max+s-1
    for i in 1:length(seq)-k+1
        strobemer = get_strobe_2_mer(view(seq,i:i+k-1), s, w_min, w_max, q; withGap = false)
        bins[as_UInt(strobemer) + 1] += 1
    end
end

# everything below here is unfinished ################################
"""
# initialzies the first initial strobemer
function initialize_first_strobe!(seq::DnaSeq, first_strobe::UInt, s::Int,
    Nt_bits::DnaBits = NUCLEOTIDE_BITS)
    for nt in view(seq, 1:s-1)
        first_strobe = (first_strobe << 2) | Nt_bits[nt]
    end
end

# just a bit where each two bits is a nt
function initialize_window!(seq::DnaSeq, window::UInt, w_min::Int, w_max::Int,
    Nt_bits::DnaBits = NUCLEOTIDE_BITS)
    for nt in view(seq, w_min:w_max-1)
        window = (window << 2) | Nt_bits[nt]
    end
end

function paired_strobe_2_mer_count(
    seq::DnaSeq,
    s::Int = 2,
    w_min::Int = 3,
    w_max::Int = 5,
    q::Int = 5;
    Nt_bits::DnaBits = NUCLEOTIDE_BITS)

    k = w_max + s - 1
    mask = unsigned(2 << (4*s - 1) - 1) # from benchmarking it was faster to create new one every time
    bins = zeros(mask + 1)
    
    first_strobe = unsigned(0); initialize_first_strobe!(seq, first_strobe, s)
    strobe_window = unsigned(0); initialize_window!(seq, strobe_window, w_min, w_max)

    for i in 1:length(seq)-k+1
        first_strobe = ((first_strobe << 2) & mask) | Nt_bits[seq[i]]
        strobemer = get_strobe_2_mer(view(seq,i:i+k-1), s, w_min, w_max, q, )
        bins[as_UInt(strobemer) + 1] += 1
    end
    return bins
end

# testing
Random.seed!(42)
vgene = dna"caggtccagctggtgcagccaggggctgagctgaggaagcctggggctttgctgaaggtctcctgcaaggcttctggatacaccttcaccagctactacatagactgggtgcgacaggcccctggacaagggcttgggtgggtgggaagaattgaccctgaagacggtggcacaaattatgcacagaagttccagggcagagtcaccttgactgcagacacgtccaccagcacagcctacgtggagctgagcagtctgagatctgaggacacggccgtgtgttactgtgtgagaga"
vgene_10_div = mutate_seq(vgene, 0.3)
rand_seq = randdnaseq(296)

seq1 = paired_strobe_2_mer_count(vgene)
seq2 = strobe_2_mer_count(vgene_10_div)
rand_seq3 = strobe_2_mer_count(rand_seq)
plot(seq1); plot!(seq2)
kmer_dist(seq1,seq2)
dist(seq1, rand_seq3)
"""