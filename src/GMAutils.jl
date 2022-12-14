"""
    getSeq(FASTA.Record)

get the dna longsequence of a fasta record, its simply calling FASTX
"""
function getSeq(seq::FASTX.FASTA.Record)
     return FASTX.FASTA.sequence(LongSequence{DNAAlphabet{4}}, seq)
 end

export getSeq

# standard dictionary based kmer frequency from scratch with some given imputs, used in GMA. 
# replacing with the proposed optimal version can save almost a minute on 4 billion bps from 5 to 4 mins
function fasterKF(k::Int64, seq::LongSubSeq{DNAAlphabet{4}}, #::LongSequence{DNAAlphabet{4}},
    KD::Dict{LongSequence{DNAAlphabet{4}}, Int64}, rv::Vector{Float64})
    k -= 1
    for i in 1:length(seq)-k
        rv[KD[view(seq, i:i+k)]] += 1 
    end
    #return rv
end

export fasterKF

#in the future i hope to be able to iterate and read at the same time through the file. That would speed it up dramatically.

"""
    genKmers(k::Int64,
             Dictionary::Bool = true;
             withN::Bool = false,
             cumulative::Bool = false)

systematic Kmer generation into dictionary：It matches each kmer to a unique index.

the withN argument can be specified to decide whether dna"N" should be included as a nucleotide in the resulting kmers.

the cumulative argument shouldnt ever be needed but it generates a kmer Dictionary with all the kmer sizes up to k.
"""
function genKmers(k::Int64, Dictionary::Bool = true; withN::Bool = false, cumulative::Bool = false) #O(4^k) no way around it im pretty sure.
    bases = [dna"A"d, dna"C"d, dna"T"d, dna"G"d]
    if withN; push!(bases,dna"N"d) end
    last = [dna""d]
    curr = LongSequence{DNAAlphabet{4}}[]
    if Dictionary == false
        for i in 1:k
            for b in bases
                for l in last
                    push!(curr,l*b)
                end
            end
            last = copy(curr) #this may be slowing it down
            curr = LongSequence{DNAAlphabet{4}}[]
        end
        return last
    else
        kmers = Dict{LongSequence{DNAAlphabet{4}}, Int64}()
        index = 0
        for i in 1:k
            for b in bases
                for l in last
                    push!(curr,l*b)
                    if !cumulative
                        if length(l*b) == k #might speed up by a few nanosecs if I keep a count then
                            index += 1
                            kmers[l*b] = index
                        end
                    else
                        index += 1
                        kmers[l*b] = index
                    end
                end
            end
            last = copy(curr)
            curr = LongSequence{DNAAlphabet{4}}[]
        end
        return kmers
    end
end

export genKmers


"""
    kmerFreq(k::Int64,
         seq::LongSequence{DNAAlphabet{4}},
         KD::Dict{LongSequence{DNAAlphabet{4}};
         returnDict::Bool = false)

Computes the kmer frequency vector of a sequence from scratch by iteration (O(n)). It is actually a slightly slower version that the main GMA does not use.

KD is the kmer Dictionary and an optional argument. If left empty, will be automatically generated WITH N's.

Note that it includes overlap(it can be tested whether it helps. but intuitively the difference shouldnt be too large).
"""
function kmerFreq(k::Int64, seq::LongSequence{DNAAlphabet{4}}, KD::Dict{LongSequence{DNAAlphabet{4}}, Int64}; returnDict::Bool = false)
    if !returnDict
        rv = fill(0.0,length(KD))
        for i in 1:length(seq)-k+1
            subseq = seq[i:i+k-1]
            if get(KD,subseq,nothing) != nothing
                rv[KD[subseq]] += 1
            end
        end
        return rv
    else
        rv = KD
        for key in rv
            rv[first(key)] = 0
        end
        for i in 1:length(seq)-k+1
            subseq = seq[i:i+k-1]
            if get(KD,subseq,nothing) != nothing
                rv[subseq] += 1
            end
        end
        return rv
    end
end

function kmerFreq(k::Int64, seq::LongSequence{DNAAlphabet{4}}; returnDict::Bool = false)
    KD = genKmers(k; withN = true)
    if !returnDict
        rv = fill(0.0,length(KD))
        for i in 1:length(seq)-k+1
            subseq = seq[i:i+k-1]
            if get(KD,subseq,nothing) != nothing
                rv[KD[subseq]] += 1
            end
        end
        return rv
    else
        rv = KD
        for key in rv
            rv[first(key)] = 0
        end
        for i in 1:length(seq)-k+1
            subseq = seq[i:i+k-1]
            if get(KD,subseq,nothing) != nothing
                rv[subseq] += 1
            end
        end
        return rv
    end
end

export kmerFreq
