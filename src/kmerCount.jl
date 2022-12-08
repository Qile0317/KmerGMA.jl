# standard dictionary based kmer frequency from scratch with some given imputs, used in GMA. 
# replacing with the proposed optimal version can save almost a minute on 4 billion bps from 5 to 4 mins
function fasterKF(k::Int64, seq::LongSubSeq{DNAAlphabet{4}}, #::LongSequence{DNAAlphabet{4}},
    KD::Dict{LongSequence{DNAAlphabet{4}}, Int64}, rv::Vector{Float64})
    k -= 1
    for i in 1:length(seq)-k
        rv[KD[view(seq, i:i+k)]] += 1 
    end
    return rv
end

export fasterKF


#heres a work in progress bit-based kmer count version that runs several orders of magnitude faster, but requires a massive re-vamp of everything
#it is based on NextGenSeqUtils
# it first requires a new type of kmer dict 

#currently theres these approaches for the bit based kmer counter that I can think of;
# 1- make it so that every time theres an N, it goes striaght to the bin and maybe divide it by k 
# 2- treat N as a nucleotide but then have to use 3 bits and needing redundant memory.
# 3- have 2 sets of dictionaries and one of them pointing valid numbers to their own new bins. this seems like it doesnt fix the issue though of hashin
 
const NUCLEOTIDE_BITS= Dict(DNA_A => unsigned(0),
                            DNA_C => unsigned(1),
                            DNA_G => unsigned(2),
                            DNA_T => unsigned(3),
                            DNA_N => unsigned(4))
       
"""
#const KmerType = Array{UInt32, 1}

#this is the approach with redundant memory but faster because the "get" function isn't needed

function kmer_count(str::LongSequence{DNAAlphabet{4}}, k::Int) #kmer, bins, and mask can be pregened later
    #initialize variables
    bins = zeros(eltype(KmerType), (2^(3*k))-1) #so unfortunately alot of bins will be redundant memory that are 0s. 
    mask = unsigned(2 << ((3*k)-1)-1)  # all ones that is the kmer length
    kmer = unsigned(0)

    #start the kmer operations
    for c in str[1:k-1]
        kmer = (kmer << 3) + NUCLEOTIDE_BITS[c]
    end
    for c in str[k:end]
        kmer = ((kmer << 3) & mask) + NUCLEOTIDE_BITS[c]
        bins[kmer + 1] += 1
    end
    return bins
end

function kmer_count(str::LongSequence{DNAAlphabet{4}}, k::Int, 
    bins::Vector, mask::UInt64) #kmer, bins, and mask can be pregened later

    kmer = unsigned(0)

    for c in str[1:k-1]
        kmer = (kmer << 3) + NUCLEOTIDE_BITS[c]
    end
    for c in str[k:end]
        kmer = ((kmer << 3) & mask) + NUCLEOTIDE_BITS[c]
        bins[kmer + 1] += 1
    end
    return bins
end

#benchmarking
KD = KmerGMA.genKmers(6,withN=true)
rv = fill(0.0,5^6)
RV = zeros(eltype(UInt32), (2^(3*6))-1)
mask = unsigned(2 << ((3*6)-1)-1)

@BenchmarkTools.benchmark kmer_count(dna"ATGCATGC",6,RV,mask) #ouch. this uses 64 times the memory. should be fine for small k though.
#128ns pm 212ns. theres alot of GC time due to how big the array is.
@BenchmarkTools.benchmark fasterKF(6,view(dna"ATGCATGC",1:8),KD,rv)
#137ns pm 23 ns...

function orig(str::LongSubSeq{DNAAlphabet{4}}, k::Int)
    # TODO: could directly encode `str` as 2-bit BioSequence
    bins = zeros(eltype(KmerType), 4^k)

    mask = unsigned(4^k - 1)  # all ones
    kmer = unsigned(0)
    for c in str[1:k-1]
        kmer = (kmer << 2) + NUCLEOTIDE_BITS[c]
    end
    for c in str[k:end]
        kmer = ((kmer << 2) & mask) + NUCLEOTIDE_BITS[c]
        bins[-~kmer] += 1
    end
    return bins
end

@benchmark orig(dna"ATGCATGC",1)
#1.3 nanoseconds... jesus christ

const Kmers = Array{UInt64, 1}

function kmer_count_new(str::LongSubSeq{DNAAlphabet{4}}, k::Int)
    #initialize variables
    bins = zeros(eltype(Kmers), 1 + 4^k) # extra bin for N's. 
    mask = unsigned(4^k - 1)  # all ones
    kmer = unsigned(0)
    scaled_N = 1 #/k
    hasN = 1

    NUCLEOTIDE_BITS= Dict(DNA_A => unsigned(0),
                            DNA_C => unsigned(1),
                            DNA_G => unsigned(2),
                            DNA_T => unsigned(3),
                            DNA_N => unsigned(4))

    # form the first kmer
    for c in str[1:k-1]
        if hasN == 1
            if c != DNA_N
                kmer = (kmer << 2) + NUCLEOTIDE_BITS[c]
            else
                bins[end] += scaled_N
                hasN = k
            end
        else
            kmer <<= 2
            hasN -= 1
            bins[end] += scaled_N
        end
    end

    # start the kmer counting process
    for c in str[k:end]
        if hasN == 1
            if c != DNA_N
                kmer = ((kmer << 2) & mask) + NUCLEOTIDE_BITS[c]
                bins[-~kmer] += 1 # -~ just adds 1, idk if its faster than +1
            else
                kmer = ((kmer << 2) & mask)
                bins[end] += scaled_N
                hasN = k 
            end
        else
            kmer <<= 2
            hasN -= 1
            bins[end] += scaled_N
        end
    end
    return bins
end

VP = open(FASTA.Reader, "C:/Users/lu_41/Desktop/Sofo Prok/VicPac32.fna")
str = KmerGMA.getSeq(first(VP))
close(VP)

@benchmark orig(view(str,1:300),6) #4.5us 
@benchmark kmer_count_new(view(str,1:300),6) #10.8 us
@BenchmarkTools.benchmark fasterKF(6,view(str,1:300),KD,rv) #15 us

#wow. adding the equating operation made it go to the GC
#if N's did not exist, then this algo would be incredibly fast.
#But having N makes it much more complicated
"""