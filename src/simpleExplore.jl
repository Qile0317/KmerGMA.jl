"""
These functions are the essential dependencies for other functions and algorithms. There are also some functions that I experimented with and thought were unessecary.
"""
#remember to have this vers of ngsu:  Pkg.add(PackageSpec(name="NextGenSeqUtils", rev="Missing-LongCharSeq-fix", url = "https://github.com/MurrellGroup/NextGenSeqUtils.jl.git"))
#all nessecary genomic testing data:
#VicPac = "C:/Users/lu_41/Desktop/Sofo Prok/VicPac32.fna"
#fVP = (first(open(FASTA.Reader,VicPac)))
#fVicPac = FASTA.sequence((first(open(FASTA.Reader,VicPac))))
#LA = "C:/Users/lu_41/Desktop/Sofo Prok/genbank_pull_lama_and_alpaca.fasta"
#V3 = "C:/Users/lu_41/Desktop/Sofo Prok/VgeneData/final/V3.fasta"
#AlpacaV = "C:/Users/lu_41/Desktop/Sofo Prok/VgeneData/AlpacaV.fasta"
#Merged25 = "C:/Users/lu_41/Desktop/Sofo Prok/VgeneData/Sequencer/Merged/BEN-25_S23_L001.vsearch_Merged.fastq"

#Functions to do with lengths and counts of sequences/assemblies, not that useful. These were a bit helpful during debugging tho
function seqsizeV(reader::FASTX.FASTA.Reader{})
    lengths = Int64[]
    for record in reader
        l = FASTX.FASTA.seqsize(record)
        push!(lengths,l)
    end
    return lengths
end

function seqsizeV(path::String)
    seqsizeV(open(FASTA.Reader, path))
end

export seqsizeV

function avgseqsize(lengths::Vector{Int64})
    return sum(lengths)/length(lengths)
end

"""
    recordCount(reader::FASTX.FASTA.Reader)

WARNING!! IT CLOSES THE READER FOR SOME REASON

counts the number of records in a FASTA.Reader object.
"""
function recordCount(reader::FASTX.FASTA.Reader)
    c = 0
    for record in reader
        c += 1
    end
    return c
end

export recordCount

#BIG PROBLEM: I THINK DOING OPEN(FASTA.READER) DOESNT WORK IN FUNCTIONS.

"""
    avgRecLen(reader;
              rnd::Bool = true)

get the rounded average length of every record in a FASTA.Reader object. (it can also be a string indicating the path)

rnd is an optional argument for whether the result should be rounded to an Interger.
"""
function avgRecLen(reader::FASTX.FASTA.Reader; rnd::Bool = true)
    Alen = 0
    d = 0
    for record in reader
        Alen += FASTX.FASTA.seqsize(record)
        d += 1
    end
    rnd ? Int64(round(Alen/d)) : Alen/d
    close(reader)
end

function avgRecLen(path::String; rnd::Bool = true)
    Alen = 0
    d = 0
    open(FASTX.FASTA.Reader, path) do reader
        for record in reader
            Alen += FASTX.FASTA.seqsize(record)
            d += 1
        end
    end
    rnd ? Int64(round(Alen/d)) : Alen/d
end

export avgRecLen

function seqMode(lengths)
    return maximum(lengths) - minimum(lengths)
end

function seqSd(reader::FASTX.FASTA.Reader)
    lengths = seqsizegths(reader)
    avg = sum(lengths)/length(lengths)
    all = []
    for length in lengths
        a = abs(length-avg)
        push!(all,a)
    end
    l = sum(all)
    return sum(all)/sqrt(length(lengths))
end

#Very important functions to do with Dictionary and vector conversions. ApproxMatch depends on these functions!

#function to flip all keys and values in a dictionary that will become useful later.
#I also should have a flipDict!() that mutates curr dict
"""
    flipDict(dict::Dict{LongSequence{DNAAlphabet{4}},Int64})

function to return a NEW dictionary that has the key value pairs flipped compared to the input.
"""
function flipDict(dict::Dict{LongSequence{DNAAlphabet{4}},Int64})
    nd = Dict{Int64,LongSequence{DNAAlphabet{4}}}()
    for key in dict
        k = first(key)
        v = last(key)
        nd[v] = k
    end
    return nd
end

export flipDict

#kmer to interger and the reverse. Much more methods are possible (e.g. removing need to dictionary) but ill implement those later.
# need to do: dict vec, vec vec, vec dict methods.
function convertSeq(Counts::Dict{LongSequence{DNAAlphabet{4}},Int64}, KD::Dict{LongSequence{DNAAlphabet{4}},Int64})
    ans = Dict{Int64,Int64}()
    for i in Counts
        i = first(i)
        ans[KD[i]] = Counts[i]
    end
    return ans
end

function convertSeq(Counts::Dict{Int64,Int64}, KmerDict::Dict{LongSequence{DNAAlphabet{4}},Int64})
    dict = flipDict(KmerDict)
    ans = Dict{LongSequence{DNAAlphabet{4}},Int64}()
    for index in Counts #remember the key is the kmer and value is the count.
        kmerInt = first(index)
        ans[dict[kmerInt]] = last(index)
    end
    return ans
end

function convertSeq(IMGTref::Dict{LongSequence{DNAAlphabet{4}},Float64}, KmerDict::Dict{LongSequence{DNAAlphabet{4}},Int64})
    ans = Dict{Int64,Float64}()
    for i in IMGTref
        i = first(i)
        ans[KmerDict[i]] = IMGTref[i]
    end
    return ans
end

#function to convert kmerdict to vector while preserving order. More methods can be added later.
function dictToVec(KmerDict::Dict{Int64, Int64})
    ans = Int64[]
    for i in 1:length(KmerDict)
        push!(ans, KmerDict[i])
    end
    return ans
end

function dictToVec(KmerDict::Dict{Int64, Float64})
    ans = Float64[]
    for i in 1:length(KmerDict)
        push!(ans, KmerDict[i])
    end
    return ans
end

#In testing, ORder IS PRESERVED! But do not use. its not the kmer dict order. cvDicVec should be used
#also might be useful to rename cvDicVec to just KFV() for kmer frequency vector.
function dictToVec(KmerDict::Dict{LongSequence{DNAAlphabet{4}}, Float64})
    ans = Float64[]
    for i in KmerDict
        push!(ans, KmerDict[first(i)])
    end
    return ans
end

#nice function to sum up 2 common functions convert to seq and dict to vec
function cvDicVec(maindict::Dict{LongSequence{DNAAlphabet{4}}, Float64}, KD::Dict{LongSequence{DNAAlphabet{4}}, Int64})
    return dictToVec(convertSeq(maindict,KD))
end

"""
    kfv(maindict::Dict{LongSequence{DNAAlphabet{4}}, Float64},
        KD::Dict{LongSequence{DNAAlphabet{4}}, Int64})

Converts a dictionary of kmers and COUNTS (maindict) to a kmer frequency VECTOR.

KD is the kmer dictionary and must be based on the the same kmer length as the count dictionary.
"""
function kfv(maindict::Dict{LongSequence{DNAAlphabet{4}}, Float64}, KD::Dict{LongSequence{DNAAlphabet{4}}, Int64})
    return cvDicVec(maindict, KD)
end

export kfv

###### misc ##########
#old function from seconddraft might be useful just counts all kmers
#alternative that does count zeroes w kmerDict. note!!! its also possible to make a version that just generates a kmerdict full of zeroes!
function recordKCount(k::Int64, seq::LongSequence{DNAAlphabet{4}},kmerDict::Dict{LongSequence{DNAAlphabet{4}},Int64}) #doesnt work for k=1 lol but that can be easily fixed later
    k = k-1
    rv = kmerDict
    for key in rv
        rv[first(key)] = 0
    end
    for i in 1:length(seq)-k
        subseq = seq[i:i+k]
        if get(rv,subseq,nothing) != nothing
            rv[subseq] += 1
        end
    end
    return rv
end

#version that does whole reader but not count zeroes. It takes way too long to run. already 5 minutes for k = 3
#there are actual good kmer counting algorithms but this is a naive implementation
function ReaderKCount(k::Int64, reader::FASTX.FASTA.Reader, kmerDict::Dict{LongSequence{DNAAlphabet{4}},Int64})
    rv = kmerDict
    for key in rv
        rv[first(key)] = 0
    end
    c = dna""
    if k > 1
        #k -= 1
        for record in reader
            seq = FASTA.sequence(record)
            Aseq = c*seq
            c = seq[end-k+2:end]
            for i in 1:length(Aseq)-k+1
                subseq = Aseq[i:i+k-1]
                if haskey(rv,subseq) #Fix by checking if its in the kmer dict.
                    rv[subseq] += 1
                end
            end
        end
    end
    return rv
end

#VP = open(FASTA.Reader,VicPac)
#@time Vp6k = ReaderKCount(6,VP,genKmers(6))


"""
    readerlens(reader::FASTX.FASTA.Reader)

function to just look at the individual record lengths of a reader object and returns a vector of lengths
"""
function readerlens(reader::FASTX.FASTA.Reader)
    lenvec = Int64[]
    for record in reader
        push!(lenvec,FASTX.FASTA.seqsize(record))
    end
    return lenvec
end

export readerlens
#needs testing

"""
    readerNTs(reader::FASTX.FASTA.Reader)

function to count the number of nts in a fasta reader object (length)
"""
function readerNTs(reader::FASTX.FASTA.Reader)
    len = 0
    for record in reader
        len += FASTX.FASTA.seqsize(record)
    end
    return len
end

function readerNTs(reader::String)
    len = 0
    open(FASTA.Reader, reader) do io
        for record in io
            len += FASTX.FASTA.seqsize(record)
        end
    end
    return len
end

export readerNTs

##looking at FASTQ sequences with weird lengths(also can be done for fasta.):
function inspectseqsize(reader::FASTX.FASTQ.Reader, lengths::UnitRange{Int64}, amount::Int64 = 10, A::Bool = false)
    c = 0
    if !A
        for record in reader
            if c <= amount
                len = FASTQ.seqsize(record)
                if len <= last(lengths) && len >= first(lengths)
                    println(record)
                    c += 1
                end
            end
        end
    else
        for record in reader
            if c <= amount
                len = FASTQ.seqsize(record)
                if len <= last(lengths) && len >= first(lengths)
                    println(">"*FASTQ.identifier(record)*" | length = "*string(len))
                    println(FASTQ.sequence(record))
                    c += 1
                end
            end
        end
    end
end

function howManySeq(reader::FASTX.FASTQ.Reader, n::Int64)
    lessthan = 0
    morethan = 0
    for record in reader
        len = FASTQ.seqsize(record)
        if len <= n
            lessthan += 1
        else
            morethan += 1
        end
    end
    return string(lessthan)*" below or equal to thr, "*string(morethan)*" above thr."
end

"""
    percentN(seq)

Simple function to count the percentage of N nucleotides in a longsequence, record, or reader object.
"""
function percentN(seq::LongSequence{DNAAlphabet{4}})
    allbp = length(seq)
    ncount = 0
    start = 1
    rg = findfirst(ExactSearchQuery(dna"N"), view(seq, start: length(seq)))
    while !isnothing(rg)
        ncount +=1
        start += last(rg)
        rg = findfirst(ExactSearchQuery(dna"N"), view(seq, start: length(seq)))
    end
    return ncount/allbp
end

function percentN(seq::FASTX.FASTA.Record)
    seq = FASTA.sequence(seq)
    percentN(seq)
end

function percentN(collection::FASTX.FASTA.Reader)
    #remember to run the length in ExactMatch first!
    allbp = length(collection)
    ncount = 0
    for seq in collection
        seq = FASTA.sequence(seq)
        start = 1
        rg = findfirst(ExactSearchQuery(dna"N"), view(seq, start: length(seq)))
        while !isnothing(rg)
            ncount +=1
            start += last(rg)
            rg = findfirst(ExactSearchQuery(dna"N"), view(seq, start: length(seq)))
        end
    end
    return ncount/allbp
    close(reader)
end

export percentN
