#script with functions needed for reference generation.
"""
    genRef(k::Int64,
       reader,
       kmerDict::Dict{LongSequence{DNAAlphabet{4}}, Int64})

Generate the reference from a FASTA.Reader object of reference sequences and returns a dictionary of the reference.

reader can be a FASTA.Reader object or a string indicating the path of the fasta file.

Is the reference generation algorithm used in the GMA.

can heavily optimized but it doesn't matter too much atm but it will probably not scale too well for very long databases
"""
function genRef(k::Int64, reader::FASTX.FASTA.Reader, kmerDict::Dict{LongSequence{DNAAlphabet{4}}, Int64}) #; returnDict::Bool = true
    len = 0
    answer = Dict{LongSequence{DNAAlphabet{4}}, Float64}()
    for key in kmerDict
        answer[LongSequence{DNAAlphabet{4}}(first(key))] = 0.0
    end
    for record in reader
        seq = getSeq(record)
        for i in 1:length(seq)-k+1
            answer[seq[i:i+k-1]] += 1
        end
        len += 1
    end
    for key in answer
        answer[first(key)] /= len
    end
    return answer
    #dont close reader!!
end

function genRef(k::Int64, path::String, kmerDict::Dict{LongSequence{DNAAlphabet{4}}, Int64})
    len = 0
    answer = Dict{LongSequence{DNAAlphabet{4}}, Float64}()
    for key in kmerDict
        answer[LongSequence{DNAAlphabet{4}}(first(key))] = 0.0
    end
    open(FASTA.Reader,path) do reader
        for record in reader
            seq = LongSequence{DNAAlphabet{4}}(FASTA.sequence(record))
            for i in 1:length(seq)-k+1
                answer[seq[i:i+k-1]] += 1
            end
            len += 1
        end
    end
    for key in answer
        answer[first(key)] /= len
    end
    return answer
end

export genRef

"""
    findthr(refseqs::Union{FASTX.FASTA.Reader, String},
            refKFV::Dict{LongSequence{DNAAlphabet{4}}, Float64},
            KD::Dict{LongSequence{DNAAlphabet{4}}, Int64};
            buff::Union{Int64,Float64} = 25)

Suggests an SED threshold, assuming most indexes do not match. Returns a float

Its extremely simple and just adds to a the SED of the first reference sequence's KFV to the actual reference KFV
"""
function findthr(refseqs::String, refKFV::Dict{LongSequence{DNAAlphabet{4}}, Float64},
    KD::Dict{LongSequence{DNAAlphabet{4}}, Int64}; buff::Union{Int64,Float64} = 25)
    answer = 0
    open(FASTX.FASTA.Reader, refseqs) do io
        seq = getSeq(first(io))
        answer += Distances.sqeuclidean(kmerFreq(
        length(first(first(KD))),seq,KD),
        kfv(refKFV,KD))
    end
    return answer + buff
end

export findthr

"""
    findRandThr(refseqs::Union{FASTX.FASTA.Reader, String},
                refKFV::Dict{LongSequence{DNAAlphabet{4}}, Float64},
                KD::Dict{LongSequence{DNAAlphabet{4}}, Int64};
                sampleSize::Int64 = 0, 
                buff::Union{Int64,Float64} = 0,
                ScaleFactor::Float64 = 0.0
                setSeed::Int64 = 1112)  

Suggests an kmer distance threshold, assuming most indexes do not match. Returns a float.

It is actually partially based on randomness as the algorithm compares to a random dna seq.
"""
function findRandThr(refseqs::String, 
    refKFV::Dict{LongSequence{DNAAlphabet{4}}, Float64},
    KD::Dict{LongSequence{DNAAlphabet{4}}, Int64};
    sampleSize::Int64 = 0, 
    buff::Union{Int64,Float64} = 0,
    ScaleFactor::Float64 = 0.0,
    setSeed::Int64 = 1112)  

    Random.seed!(setSeed) 
    k = length(first(first(KD))) #get kmer size from kmer dict
    if ScaleFactor == 0.0; ScaleFactor = 1/(2*k) end #set scale factor
    if sampleSize == 0; sampleSize = 300 end #set samplesize

    io = open(FASTX.FASTA.Reader, refseqs)
    seq = first(io) #get sequence of the first record of reference
    close(io)

    seqlen = FASTA.seqsize(seq)
    seq = getSeq(seq)

    #find the distance of the first reference to the KFV, should be small. 
    exact = ScaleFactor * Distances.sqeuclidean(
        kmerFreq(k,seq,KD),kfv(refKFV,KD)) 

    upp = 0
    #Random.seed!(setSeed) #set seed for random sequence generation. ill add in next updates 

    for i in 1:sampleSize #run n random sequences 
        upp += ScaleFactor * Distances.sqeuclidean(
            kmerFreq(k,BioSequences.randdnaseq(seqlen),
            KD),kfv(refKFV,KD)) #check against a random sequence
    end
    upp = upp * ScaleFactor / 300
    return exact + ((upp-exact)/2) + buff #returns halfway point. It should be mathematically investigated in the future whether the halfway point is optimal
end

export findRandThr
