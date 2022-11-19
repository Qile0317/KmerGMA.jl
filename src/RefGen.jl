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

Suggests an SED threshold, assuming most indexes do not match. Returns an interger

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
#findthr should be able to just take KFV instead of requireing KFD. its an ez fix tho ill do it later
function findthr(refseqs::FASTX.FASTA.Reader,
    refKFV::Dict{LongSequence{DNAAlphabet{4}}, Float64},
    KD::Dict{LongSequence{DNAAlphabet{4}}, Int64};
    buff::Union{Int64,Float64} = 25)

    seq = getSeq(first(refseqs))
    answer = Distances.sqeuclidean(kmerFreq(length(first(first(KD))),
    seq,KD), kfv(refKFV,KD))
    return answer + buff
end

export findthr

"""
    findavgthr(refseqs::Union{FASTX.FASTA.Reader, String},
               refKFV::Dict{LongSequence{DNAAlphabet{4}}, Float64},
               KD::Dict{LongSequence{DNAAlphabet{4}}, Int64};
               buff::Union{Int64,Float64} = 25)

version of findthr that scans through the entire reference to see the average SED for the most accurate avg SED but probably doesnt make much of a difference.
"""
function findavgthr(refseqs::Union{FASTX.FASTA.Reader, String}, refKFV::Dict{LongSequence{DNAAlphabet{4}}, Float64},
    KD::Dict{LongSequence{DNAAlphabet{4}}, Int64}; buff::Union{Int64,Float64} = 25)
    if typeof(refseqs) == String
        refseqs = open(FASTA.Reader, refseqs)
    end
    SED = []
    KFV = kfv(refKFV,KD)
    k = length(first(first(KD)))
    for record in refseqs
        seq = FASTA.sequence(record)
        push!(SED, Distances.sqeuclidean(kmerFreq(k,seq,KD),KFV))
    end
    return (sum(SED)/recordCount(refseqs)) + buff
end

#export findavgthr
