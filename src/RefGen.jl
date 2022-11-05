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
    len = recordCount(reader) # this screws the function
    answer = Dict{LongSequence{DNAAlphabet{4}}, Float64}()
    for key in kmerDict
        answer[first(key)] = 0.0
    end
    for record in reader
        seq = FASTA.sequence(record)
        for i in 1:length(seq)-k+1
            answer[seq[i:i+k-1]] += 1
        end
    end
    for key in answer
        answer[first(key)] /= len
    end
    return answer
end

function genRef(k::Int64, path::String, kmerDict::Dict{LongSequence{DNAAlphabet{4}}, Int64})
    reader = open(FASTA.Reader,path)
    len = recordCount(reader)
    answer = Dict{LongSequence{DNAAlphabet{4}}, Float64}()
    for key in kmerDict
        answer[first(key)] = 0.0
    end
    reader = open(FASTA.Reader,path)
    for record in reader
        seq = FASTA.sequence(record)
        for i in 1:length(seq)-k+1
            answer[seq[i:i+k-1]] += 1
        end
    end
    for key in answer
        answer[first(key)] /= len
    end
    return answer
    close(reader)
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
function findthr(refseqs::Union{FASTX.FASTA.Reader, String}, refKFV::Dict{LongSequence{DNAAlphabet{4}}, Float64},
    KD::Dict{LongSequence{DNAAlphabet{4}}, Int64}; buff::Union{Int64,Float64} = 25)
    if typeof(refseqs) == String
        refseqs = open(FASTA.Reader, refseqs)
    end
    seq = FASTA.sequence(first(refseqs))
    return Distances.sqeuclidean(kmerFreq(length(first(first(KD))),seq,KD),
    kfv(refKFV,KD)) + buff
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

export findavgthr