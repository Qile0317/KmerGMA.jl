"""
TO DO/TO ASK LIST:

        -> big issue: the function doesnt account for overlap between records in a reader. can be implemented in the future tho.
        -> make a threshold predictor by taking averages and going 1-2 SD down. although for k=6, usually around SED = 150 - 190 works ok?
        -> There can be a currnt minimum function that returns it so there doesnt have to be so many repeated code chunks of the first 2 operations
            -> alternative approach is to take any sequence in reference and compare it to the reference KFV and take that as average.
        -> the buffer thing needs to be fixed if it exceeds
        -> implement progress bar.
        -> make plotting function for where the matches were
        -> keep track of global position with a count variable.

        -> potentially changing the int type to unsigned might speed up the code

        -> next step in pipeline is to take the buffers, align, and trim off excess ends and then unalign
        -> then webBLAST

        -> theres also a potentially completely unessecary step to cluster similar reference sequences together.

        -> N stuff:
            -> N THRESHOLD!!! (interger and truncated percentage!) maybe. probably a bad idea actually.
            -> make a function tha detects start and end points of Ns in the genome and plot it. Just to see if it is covering ay important regions.
        -> implement bitwise operations & euclidian dist to speed up code (probably ditching)
        -> write a master function that does this from beginning to end in one super long and slow step.
        -> FOR PAPER: CHECK DIFFERENT K VALUES AND EFFECTIVENESS/IF IT EXISTS. i THINK PROBABLY 4-6 IS THE OPTIMAL RATIO BUT USE IT TO SEE SPARSITY AND MATCHES ACCURACY.

        ->and ofc compare to state of the art HMM models. try to one up them on just speed

        -> in final funcition have an argument omit::Bool to check whether or not to omit exact matches.

        -> !!!! COMPARE THIS METHOD WITH EXACTMATCH.JL AT Euclidian distance 0.

        -> implement phylogeny stuff like ben.
            -> on paper, put phylogeny paired with eucdist

        -> Implement golden formula later.

        -> generate global kmer dict consts at the start. can actually utilize the genKmers()
           function cut off the same kmer.

        ->!!!!!!!!!!!!!
            TEST EACH ALPACA SEQ AGAINST AVERAGE IMGT REF!!! IT MIGHT BE LIKE 18 ISH!!!!

        -> make it possible to generate reference with an alignment as well
        -> For the future, in the eucvector you should be able to see how far it is from the minimum
        -> make jupyter notebook pipeline that is eaier to showcase, If I have time, showcase in real time in the presentation

In Progress:
        -> the dictionary order as of rn is Sequence => interger. Should it be reversed for the freq table?
        -> use IMGT fasta file fownloads to make the frequency table. calculate the kmer distribution for each sequence and average the kmer distributions. Generate the data and store!!!
        -> implement rev in some kmer count functions in order to reverse dictionary key order.
        -> Kmer spectra for the future.
"""
function toDoList()
    println("type ?toDoList in REPL")
end

#using Pkg
#Pkg.add(PackageSpec(name="NextGenSeqUtils", rev="Missing-LongCharSeq-fix", url = "https://github.com/MurrellGroup/NextGenSeqUtils.jl.git"))

#using BioSequences, FASTX, Plots, ProgressMeter, Distances, PkgTemplates, Pkg, DataFrames, Distances, DelimitedFiles#, Distributed


#implementation of functions needed for the first step of the optimized method
#frequency vector with no actual kmer information. Should only be ran once in the master function
#order matters. goes from left to right in the order of the kmer.
#it is O(n) assuming hashing is almost always O(1).
 #reference generation depends on the function kmerFreq!

#faster kmer frequency with some given imputs.
 function fasterKF(k::Int64, seq::LongSequence{DNAAlphabet{4}}, KD::Dict{LongSequence{DNAAlphabet{4}}, Int64}, rv::Vector{Float64})
     for i in 1:length(seq)-k+1
         subseq = seq[i:i+k-1]
         if get(KD,subseq,nothing) != nothing
             rv[KD[subseq]] += 1
         end
     end
     return rv
 end

#this version just treates N as another nucleotide. it actually seems to work well as Ns incur a high sqeuclidean dist.
function queryMatch(k::Int64, record::FASTX.FASTA.Record, IMGTref::Dict{LongSequence{DNAAlphabet{4}}, Float64}, kmerDict::Dict{LongSequence{DNAAlphabet{4}}, Int64}, windowsize::Int64 = 0)
    if windowsize == 0
        windowsize = avgRecLen(IMGTref)
    end
    seq = FASTA.sequence(record)
    IMGTrefVec = cvDicVec(IMGTref,kmerDict)
    lkd = length(kmerDict)
    rv = fill(0.0,lkd)
    curr = kmerFreq(k,seq[1:windowsize],kmerDict)

    currSqrEuc = Distances.sqeuclidean(IMGTrefVec,curr) #CURR IS A DICTUINARY
    EucVec = Float64[currSqrEuc]

    #iteration with sliding window
    for i in 2:(length(seq)-windowsize)
        #first operation

        #zeroth kmer
        zerokInt = kmerDict[seq[i-1:i+k-2]] #if i just shift everything left it'd save the subtraction time right?
        MinusOld = (IMGTrefVec[zerokInt]-curr[zerokInt])^2
        ZeroOld = curr[zerokInt]
        curr[zerokInt] -= 1
        MinusNew = (IMGTrefVec[zerokInt]-curr[zerokInt])^2

        #last kmer
        KendInt = kmerDict[seq[i+windowsize-k+1:i+windowsize]]
        EndOld = curr[KendInt]
        curr[KendInt] += 1
        PlusNew = (IMGTrefVec[KendInt]-curr[KendInt])^2
        PlusOld = (IMGTrefVec[KendInt]-EndOld)^2

        #second operation. no third
        currSqrEuc += PlusNew + MinusNew - PlusOld - MinusOld
        push!(EucVec, currSqrEuc)
    end
    return EucVec
end

export queryMatch

#incorporates faster kmer match
function queryMatchN(k::Int64, record::FASTX.FASTA.Record,
    IMGTref::Dict{LongSequence{DNAAlphabet{4}}, Float64},
    kmerDict::Dict{LongSequence{DNAAlphabet{4}}, Int64}, windowsize::Int64,
    lengthkd::Int64, rv::Vector{Float64})

    seq = FASTA.sequence(record)
    IMGTrefVec = cvDicVec(IMGTref,kmerDict)
    curr = fasterKF(k,seq[1:windowsize],kmerDict,rv)

    currSqrEuc = Distances.sqeuclidean(IMGTrefVec,curr) #CURR IS A DUCTUINARY
    EucVec = Float64[currSqrEuc]

    #iteration with sliding window
    for i in 2:(length(seq)-windowsize)
        #first operation

        #zeroth kmer
        zerokInt = kmerDict[seq[i-1:i+k-2]] #unessecary lol
        MinusOld = (IMGTrefVec[zerokInt]-curr[zerokInt])^2
        ZeroOld = curr[zerokInt]
        curr[zerokInt] -= 1
        MinusNew = (IMGTrefVec[zerokInt]-curr[zerokInt])^2

        #last kmer
        KendInt = kmerDict[seq[i+windowsize-k+1:i+windowsize]]
        EndOld = curr[KendInt]
        curr[KendInt] += 1
        PlusNew = (IMGTrefVec[KendInt]-curr[KendInt])^2
        PlusOld = (IMGTrefVec[KendInt]-EndOld)^2

        #second operation. no third
        currSqrEuc += PlusNew + MinusNew - PlusOld - MinusOld
        push!(EucVec, currSqrEuc)
    end
    return EucVec
end

#version that prints to repl by finding minima
function printQueryMatchN(k::Int64, record::FASTX.FASTA.Record, IMGTref::Dict{LongSequence{DNAAlphabet{4}}, Float64}, kmerDict::Dict{LongSequence{DNAAlphabet{4}}, Int64} = nothing,  thr::Union{Float64, Int64} = 180.0,
    windowsize::Int64 = nothing, buff::Int64 = 50)
    thrBuf = " | thr = "*string(thr)*" | buffer = "*string(buff)
    if isnothing(windowsize)
        windowsize = avgRecLen(IMGTref)
    end

    if isnothing(kmerDict)
        kmerDict = genKmers(k,withN=true)
    end

    seq = FASTA.sequence(record)
    IMGTrefVec = cvDicVec(IMGTref,kmerDict)
    curr = kmerFreq(k,seq[1:windowsize],kmerDict,false)
    currSqrEuc = Distances.sqeuclidean(IMGTrefVec,curr)

    CMI = 1
    stop = true
    currminim = currSqrEuc

    for i in 2:(length(seq)-windowsize)
        #first operation
        #zeroth kmer
        zerokInt = kmerDict[seq[i-1:i+k-2]] #unessecary lol
        MinusOld = (IMGTrefVec[zerokInt]-curr[zerokInt])^2
        ZeroOld = curr[zerokInt]
        curr[zerokInt] -= 1
        MinusNew = (IMGTrefVec[zerokInt]-curr[zerokInt])^2

        #last kmer
        KendInt = kmerDict[seq[i+windowsize-k+1:i+windowsize]]
        EndOld = curr[KendInt]
        curr[KendInt] += 1
        PlusNew = (IMGTrefVec[KendInt]-curr[KendInt])^2
        PlusOld = (IMGTrefVec[KendInt]-EndOld)^2

        #second operation.
        currSqrEuc += PlusNew + MinusNew - PlusOld - MinusOld

        #third operation
        if currSqrEuc < thr
            if currSqrEuc < currminim
                currminim = currSqrEuc
                CMI = i
                stop = false
            end
        else
            if stop == false
                println(">"*FASTA.identifier(record)*"| SED = "*string(currminim)*" | MatchPos = "*string(CMI)*":"*string(CMI+windowsize)*thrBuf) #i can make it CMI - CMI+windowsize
                if i-buff > 0 #lol doesnt account for if buffer exceeds end.
                    println(seq[i-buff:i+windowsize-1+buff])
                else #this probably also adds time, ill make a slightly faster version later.
                    println(seq[i:i+windowsize-1+buff])
                end
                currminim = currSqrEuc
                stop = true
            end
        end
    end
end

#now for the entire reader by printing to repl
function printQueryMatchN(k::Int64, reader::FASTX.FASTA.Reader, IMGTref::Dict{LongSequence{DNAAlphabet{4}}, Float64}, kmerDict::Dict{LongSequence{DNAAlphabet{4}}, Int64} = nothing,  thr::Union{Float64, Int64} = 180.0,
    windowsize::Int64 = nothing, buff::Int64 = 50)
    #lr = recordCount(reader)
    #c = 0
    for record in reader
        #c += 1
        printQueryMatchN(k, record, IMGTref, kmerDict, thr, windowsize, buff)
    end
end

#printQueryMatchN(6, open(FASTA.Reader,LA), V3NRef, genKmers(6,withN=true), 180, 290, 30)

"""
    writeQueryMatch(k::Int64,
                    record::FASTX.FASTA.Record,
                    IMGTref::Dict{LongSequence{DNAAlphabet{4}}, Float64},
                    kmerDict::Dict{LongSequence{DNAAlphabet{4}}, Int64},
                    thr::Union{Float64, Int64} = 180.0,
                    windowsize::Int64 = 0,
                    buff::Int64 = 50,
                    path::String = nothing)

A slighty slower query match that writes to a fasta file for a SINGLE FASTA.Record. An optimized version is used in the actual GMA.
"""
function writeQueryMatch(k::Int64, record::FASTX.FASTA.Record, IMGTref::Dict{LongSequence{DNAAlphabet{4}}, Float64}, kmerDict::Dict{LongSequence{DNAAlphabet{4}}, Int64},  thr::Union{Float64, Int64} = 180.0,
    windowsize::Int64 = 0, buff::Int64 = 50, path::String = nothing)
    if buff != 0
        thrBuf = " | thr = "*string(thr)*" | buffer = "*string(buff)
    else
        thrBuf = " | thr = "*string(thr)
    end
    if windowsize == 0
        windowsize = avgRecLen(IMGTref)
    end
    seq = FASTA.sequence(record)
    IMGTrefVec = cvDicVec(IMGTref,kmerDict)
    curr = kmerFreq(k,seq[1:windowsize],kmerDict,false)
    currSqrEuc = Distances.sqeuclidean(IMGTrefVec,curr)

    CMI = 1
    stop = true
    currminim = currSqrEuc

    for i in 2:(length(seq)-windowsize)
        #first operation
        #zeroth kmer
        zerokInt = kmerDict[seq[i-1:i+k-2]] #unessecary lol
        MinusOld = (IMGTrefVec[zerokInt]-curr[zerokInt])^2 #i think instead of squaring, doing a bitshift is faster.
        ZeroOld = curr[zerokInt]
        curr[zerokInt] -= 1
        MinusNew = (IMGTrefVec[zerokInt]-curr[zerokInt])^2

        #last kmer
        KendInt = kmerDict[seq[i+windowsize-k+1:i+windowsize]]
        EndOld = curr[KendInt]
        curr[KendInt] += 1
        PlusNew = (IMGTrefVec[KendInt]-curr[KendInt])^2
        PlusOld = (IMGTrefVec[KendInt]-EndOld)^2

        #second operation.
        currSqrEuc += PlusNew + MinusNew - PlusOld - MinusOld

        #third operation
        if currSqrEuc < thr
            if currSqrEuc < currminim
                currminim = currSqrEuc
                CMI = i
                stop = false
            end
        else
            if stop == false #in future account for if buffer exceeds end or front
                rec = FASTA.Record(FASTA.identifier(record), "| SED = "*string(currminim)[1:5]*" | Pos = "*string(CMI)*":"*string(CMI+windowsize)*thrBuf, seq[i-buff:i+windowsize-1+buff])
                w = FASTA.Writer(open("VicPacScan/vicpacscan.fasta", "a"), width = 95)
                write(w, rec)
                close(w)

                currminim = currSqrEuc
                stop = true
            end
        end
    end
end

export writeQueryMatch

#faster one, essential for GMA
function writeQueryMatchN(k::Int64, record::FASTX.FASTA.Record,
    IMGTrefVec::Vector{Float64},
    kmerDict::Dict{LongSequence{DNAAlphabet{4}}, Int64},
    thr::Union{Float64, Int64} = 180.0,
    windowsize::Int64 = 0, buff::Int64 = 50, path::String = nothing;
    rv::Vector{Float64}, thrbuff::String)

    #initial operations
    seq = FASTA.sequence(record)
    curr = fasterKF(k,seq[1:windowsize],kmerDict,rv)
    currSqrEuc = Distances.sqeuclidean(IMGTrefVec,curr)

    CMI = 1
    stop = true
    currminim = currSqrEuc

    for i in 2:(length(seq)-windowsize)
        #first operation
        #zeroth kmer
        zerokInt = kmerDict[seq[i-1:i+k-2]] #unessecary lol
        MinusOld = (IMGTrefVec[zerokInt]-curr[zerokInt])^2
        curr[zerokInt] -= 1
        MinusNew = (IMGTrefVec[zerokInt]-curr[zerokInt])^2

        #last kmer
        KendInt = kmerDict[seq[i+windowsize-k+1:i+windowsize]]
        PlusOld = (IMGTrefVec[KendInt]-curr[KendInt])^2
        curr[KendInt] += 1
        PlusNew = (IMGTrefVec[KendInt]-curr[KendInt])^2

        #second operation.
        currSqrEuc += PlusNew + MinusNew - PlusOld - MinusOld

        #third operation
        if currSqrEuc < thr
            if currSqrEuc < currminim
                currminim = currSqrEuc
                CMI = i
                stop = false
            end
        else
            if stop == false #in future account for if buffer exceeds end or front
                rec = FASTA.Record(FASTA.identifier(record), "| SED = "*string(currminim)[1:5]*" | Pos = "*string(CMI)*":"*string(CMI+windowsize)*thrbuff, seq[i-buff:i+windowsize-1+buff])
                w = FASTA.Writer(open("VicPacScan/vicpacscan.fasta", "a"), width = 95) #v important that its "a"
                write(w, rec)
                close(w)

                currminim = currSqrEuc
                stop = true
            end
        end
    end
end

#testing. WORKS!!!
#d = KmerGMA.genKmers(6,withN=true)
#r = KmerGMA.genRef(6,"C:/Users/lu_41/Desktop/Sofo Prok/VgeneData/AlpacaV.fasta",d)
#ref = KmerGMA.kfv(r,d)
#rV = fill(0.0,5^6)
#red = open(FASTA.Reader,KmerGMA.LA)
#@time KmerGMA.writeQueryMatchN(6,first(red),ref,d,200.0,289,50,"VicPacScan/vicpacscan.fasta"; rv = rV, thrbuff = "test")

#using FlameGraphs, ProfileView, Profile
#Profile.clear(); @profile kmerGMA.writeQueryMatchN(6,LA,ref,d,190.0,289,50,"VicPacScan/vicpacscan.fasta"; rv = rV, thrbuff = "test")
#g = flamegraph()
#ProfileView.view()
"""
    writeQueryMatch(k::Int64,
                    reader::FASTX.FASTA.Reader,
                    IMGTref::Dict{LongSequence{DNAAlphabet{4}}, Float64},
                    path::String;
                    kmerDict::Union{Dict{LongSequence{DNAAlphabet{4}}, Int64},Nothing} = nothing,
                    thr::Union{Float64, Int64} = 180.0,
                    windowsize::Int64 = 0,
                    buff::Int64 = 50)

The main implementation of the genome mining algorithm as described by the paper.

- unfinished docs
"""
function writeQueryMatch(k::Int64, reader::FASTX.FASTA.Reader,
    IMGTref::Dict{LongSequence{DNAAlphabet{4}}, Float64}, path::String,
    kmerDict::Union{Dict{LongSequence{DNAAlphabet{4}}, Int64},Nothing} = nothing,
    thr::Union{Float64, Int64} = 180.0, windowsize::Int64 = 0, buff::Int64 = 50)

    #initializing nessecary arguments
    RV = fill(0.0,length(kmerDict))

    if windowsize == 0
        windowsize = avgRecLen(IMGTref)
    end

    IMGTref = cvDicVec(IMGTref,kmerDict) #this is a bit dumb, its actually possible to make it generate a vector in the first place.

    if buff != 0
        thrBuf = " | thr = "*string(thr)*" | buffer = "*string(buff)
    else
        thrBuf = " | thr = "*string(thr)
    end

    if isnothing(kmerDict)
        kmerDict = genKmers(k,withN=true)
    end

    #Genome scanning by writing to a file
    for record in reader
        writeQueryMatchN(k, record, IMGTref, kmerDict, thr, windowsize, buff, path;
        rv = RV, thrbuff = thrBuf)
    end
    close(reader)
end

export writeQueryMatch

#@time writeQueryMatch(6,open(FASTA.Reader,LA),V3NRef,sixMerNDict,170.0,289,50,"VicPacScan/vicpacscan.fasta")
# 0.158826 seconds (1.94 M allocations: 97.667 MiB)

#@time kmerGMA.writeQueryMatch(6,open(FASTA.Reader,LA),V3NRef,sixMerNDict,170.0,289,50,"VicPacScan/vicpacscan.fasta")


#scanning all of vicpac. i should put a progress bar.
#@time writeQueryMatch(6,open(FASTA.Reader,VicPac),V3NRef,genKmers(6,withN=true),350,289,50,"VicPacScan/vicpacscan.fasta")
#@time didnt work but this run started at around 7:06 and ended at 7:23 so it took 17 minutes ish.
