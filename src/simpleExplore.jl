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
    if withN == false
        bases = [dna"A"d, dna"C"d, dna"T"d, dna"G"d]
    else
        bases = [dna"A"d, dna"C"d, dna"T"d, dna"G"d, dna"N"d]
    end
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

export genKmers, sixMerNDict

#Functions to do with lengths and counts of sequences/assemblies, not that useful
function seqLenV(reader::FASTX.FASTA.Reader{})
    lengths = Int64[]
    for record in reader
        l = FASTX.FASTA.seqlen(record)
        push!(lengths,l)
    end
    return lengths
end

function seqLenV(path::String)
    seqLenV(open(FASTA.Reader, path))
end

function avgSeqLen(lengths::Vector{Int64})
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
    avgRecLen(reader::FASTX.FASTA.Reader, rnd::Bool = false)

get the rounded average length of every record in a FASTA.Reader object. (it can also be a string indicating the path)

rnd is an optional argument for whether the result should be rounded to an Interger.
"""
function avgRecLen(reader::FASTX.FASTA.Reader, rnd::Bool = true)
    Alen = 0
    for record in reader
        Alen += FASTX.FASTA.seqlen(record)
    end
    reader = open(FASTA.Reader, path)
    d = recordCount(reader)
    if rnd
        return Int64(round(Alen/d))
    else
        return Alen/d
    end
end

function avgRecLen(path::String, rnd::Bool = true)
    reader = open(FASTA.Reader, path)
    Alen = 0
    for record in reader
        Alen += FASTX.FASTA.seqlen(record)
    end
    if rnd == true
        reader = open(FASTA.Reader, path)
        d = recordCount(reader)
        a = Alen/d
        return Int64(round(a))
    else
        reader = open(FASTA.Reader, path)
        d = recordCount(reader)
        return Alen/d
    end
end

export avgRecLen

function seqMode(lengths)
    return maximum(lengths) - minimum(lengths)
end

function seqSd(reader::FASTX.FASTA.Reader)
    lengths = SeqLengths(reader)
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

"""
    refPlot(k::Int64,
    referece::Vector{Float64};
    nonzero::Bool = false)

A not so useful function to plot a kmer dictionary with some labels. One can easily customize the plot with their own plotting functions.
"""
function refPlot(k::Int64, referece::Vector{Float64}; nonzero::Bool = false)
    if nonzero == false
        plot(collect(1:1:length(referece)),referece,label=nothing)
        xlabel!(string("All Unique ",k,"-mers"))
        ylabel!("Average counts")
        title!(string("Average ",k,"-mer distribution in the reference"))
    elseif nonzero == true
        v = Int64[]
        for i in referece
            if i == 0.0
                push!(v,0)
            else
                push!(v,1)
            end
        end
        scatter(v,label=nothing)
        xlabel!(string("All Unique ",k,"-mers"))
        ylabel!("occurence")
        title!(string("all ",k,"-mer occurences in the average of the reference"))

    end
end

function refPlot(referece::Dict{LongSequence{DNAAlphabet{4}}, Float64}; nonzero::Bool = false, k::Int64 = 0)
    if k==0
        k = length(first(first(V3Ref)))
    end
    referece = dictToVec(referece)
    if nonzero == false
        plot(collect(1:1:length(referece)),referece,label=nothing)
        xlabel!(string("All Unique ",k,"-mers"))
        ylabel!("Average counts")
        title!(string("Average ",k,"-mer distribution in the reference"))
    elseif nonzero == true
        v = Int64[]
        for i in referece
            if i == 0.0
                push!(v,0)
            else
                push!(v,1)
            end
        end
        scatter(v,label=nothing)
        xlabel!(string("All Unique ",k,"-mers"))
        ylabel!("occurence")
        title!(string("all ",k,"-mer occurences in the average of the reference"))
    end
end

export refPlot

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

##Plotting results (I also need to do a different type of plot for eucledian distances.)
function freqTable(k::Int64, FreqDict::Dict{Int64,Int64}, Sort::Bool = true) #I need to make the title where you can change it according to k. Ik theres an easy way to do it but i havent looked it up yet.
    if Sort == true
        bar(sort(collect(keys(FreqDict))), collect(values(FreqDict)), orientation=:vertical, label = nothing)
    elseif Sort == false
        bar(collect(keys(Freqs)), collect(values(FreqDict)), orientation=:vertical, label = nothing)
    end
    xlabel!("All Unique kmers")
    ylabel!("Kmer Counts")
    title!(string(k,"-mer Counts in sequence"))
end

#adding method for vector.
function freqTable(FreqVec::Vector{Int64})
    k = Int64(round(sqrt(sqrt(length(FreqVec)))))
    bar((FreqVec), orientation=:vertical, label = nothing)
    xlabel!("All Unique kmers")
    ylabel!("Kmer Counts")
    title!(string(k,"-mer Counts in sequence"))
end

#I also need to plot the euclidian distances and some otherstuff
# I should alos implement an optional "threshold option to just plot a line"
function eucTable(eucVec::Vector{Float64}, thr::Union{Float64,Int64} = 0)
    e = plot(eucVec, label = nothing,linewidth=0.3)
    if thr != 0
        plot(e, fill(minimum(eucVec)+thr,length(eucVec)),label = "Threshold")
    else
        plot(eucVec, label = nothing,linewidth=0.3)
    end
    xlabel!("Position in genome(bp)")
    ylabel!("Euclidian distance to reference")
end

#maybe kmer spectra that fills a vector. KAT.jl has cooler ways to do this.
#problem: for something like VicPac theres too many different values. Perhaps making a histogram with each value or making some CDF is better.
function kmerSpectra(k::Int64, KCDict::Dict{LongSequence{DNAAlphabet{4}},Int64}, rm::Bool = false)
    ans = fill(0,(2^((4*k-1)))-(2^((2*k)-1)))
    for key in KCDict
        ans[last(key)+1] += 1
    end
    if rm
        Ransvec = reverse(ans)
        stop = false
        for i in Ransvec
            if stop == false
                if i == 0
                    pop!(ans)
                else
                    stop = true #there should be a better way
                end
            else
                break
            end
        end
    end
    Plots.scatter(ans,label=nothing)
    xlabel!(string(k)*"-mer frequency + 1")
    ylabel!("count")
end

function oldKmerSpectra(k::Int64, KCDict::Dict{LongSequence{DNAAlphabet{4}},Int64},typ::String="bar")
    ans = Dict{Int64,Int64}()
    @showprogress "initializing frequency array..." for key in KCDict
        if haskey(ans,last(key))
            ans[last(key)] += 1
        else
            ans[last(key)] = 1
        end
    end
    ansvec = fill(0,(2^((4*k-1)))-(2^((2*k)-1))) #this is to minimize the amount of terms to put in the vector
    for key in ans
        ansvec[first(key)+1]=last(key)
    end
    Ransvec = reverse(ansvec)
    stop = false
    for i in Ransvec
        if stop == false
            if i == 0
                pop!(ansvec)
            else
                stop = true #there should be a better way
            end
        else
            break
        end
    end
    if typ == "bar"
        Plots.bar(ansvec,label=nothing)
    else
        Plots.scatter(ansvec,label=nothing)
    end
    xlabel!(string(k)*"-mer frequency + 1")
    ylabel!("count")
end

#oldKmerSpectra(6,fc)
#oldKmerSpectra(6,Vp6k)

"""
    readerlens(reader::FASTX.FASTA.Reader; plt::Bool = true)

function to just look at the individual record lengths of a reader object and optionally plot them with the plt argument.
"""
function readerlens(reader::FASTX.FASTA.Reader; plt::Bool = true)
    lenvec = Int64[]
    for record in reader
        push!(lenvec,FASTX.FASTA.seqlen(record))
    end
    if !plt
        return lenvec
    else
        Plots.bar(lenvec, title="lengths of individual records in the FASTA file", label=nothing)
        xlabel!("sequence index")
        ylabel!("sequence length")
    end
end

export readerlens

"""
    readerNTs(reader::FASTX.FASTA.Reader)

function to count the number of nts in a fasta reader object (length)
"""
function readerNTs(reader::FASTX.FASTA.Reader)
    len = 0
    for record in reader
        len += FASTX.FASTA.seqlen(record)
    end
    return len
end

export readerNTs

##looking at FASTQ sequences with weird lengths(also can be done for fasta.):
function inspectSeqLen(reader::FASTX.FASTQ.Reader, lengths::UnitRange{Int64}, amount::Int64 = 10, A::Bool = false)
    c = 0
    if !A
        for record in reader
            if c <= amount
                len = FASTQ.seqlen(record)
                if len <= last(lengths) && len >= first(lengths)
                    println(record)
                    c += 1
                end
            end
        end
    else
        for record in reader
            if c <= amount
                len = FASTQ.seqlen(record)
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
        len = FASTQ.seqlen(record)
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
