#using Pkg
#Pkg.add(PackageSpec(name="NextGenSeqUtils", rev="Missing-LongCharSeq-fix", url = "https://github.com/MurrellGroup/NextGenSeqUtils.jl.git"))

#this version just treates N as another nucleotide. it actually seems to work well as Ns incur a high sqeuclidean dist. (actually i just made it so that it doesnt. so lets see what happens)
"""
    gma(; k::Int64,
          record::FASTX.FASTA.Record,
          refVec::Vector{Float64},
          windowsize::Int64,
          kmerDict::Dict{LongSequence{DNAAlphabet{4}}, Int64},
          thr::Float64,
          buff::Int64,
          rv::Vector{Float64},
          thrbuff::String,
          mode::String = "write",
          path::String = "noPath"
          resultVec::Vector{FASTA.Record} = FASTA.Record[],
          ScaleFactor::Float64 = 1.0)

for devs only.

this is the main genome mining algorithm for a single record. either prints, or edit a vector or file in place

the ```mode``` argument has 3 valid imputs: "write", "return", "print"
"""
function gma(;
    k::Int64,
    record::FASTX.FASTA.Record,
    refVec::Vector{Float64},
    windowsize::Int64,
    kmerDict::Dict{LongSequence{DNAAlphabet{4}}, Int64},
    thr::Float64,
    buff::Int64,
    rv::Vector{Float64},
    thrbuff::String,
    mode::String = "write",  #"return", "print"
    path::String = "noPath",
    resultVec::Vector{FASTA.Record} = FASTA.Record[],
    ScaleFactor::Float64 = 1.0) #1/2k

    sl = FASTX.FASTA.seqsize(record)
    if sl >= windowsize
        #initial operations for the first window 
        seq = getSeq(record) #getSeq is kinda slow. Ideally I wanna work with subseqs eventually. Even better is if I can read and edit at the same time
        curr = fasterKF(k,view(seq,1:windowsize),kmerDict,rv) #Theres an even faster way with unsigned ints and bits. I wonder if making so many temp arrays is a good idea
        currSqrEuc = Distances.sqeuclidean(refVec,curr)

        #initializing certain variables
        CMI = 2
        stop = true
        currminim = currSqrEuc

        for i in 1:(sl-windowsize) 
            #first & second operation
            #zeroth kmer
            @views zerokInt = kmerDict[view(seq,i:i+k-1)] #might be slightly faster to predefine
            @views currSqrEuc -= (refVec[zerokInt]-curr[zerokInt])^2 # a_old
            curr[zerokInt] -= 1
            @views currSqrEuc += (refVec[zerokInt]-curr[zerokInt])^2 # a_new

            #last kmer
            @views KendInt = kmerDict[view(seq,i+windowsize-k:i+windowsize-1)]
            @views currSqrEuc -= (refVec[KendInt]-curr[KendInt])^2 # b_old
            curr[KendInt] += 1
            @views currSqrEuc += (refVec[KendInt]-curr[KendInt])^2 # b_new

            #convert to kmer Distance 
            kmerDist = currSqrEuc * ScaleFactor 

            #third operation - the block below !stop's speed is irrelevant as its rare
            if kmerDist < thr
                if kmerDist < currminim
                    currminim = kmerDist
                    CMI = i
                    stop = false
                end
            elseif !stop #in future account for if buffer exceeds end or front
                #create the record of the match
                rec = FASTA.Record(String(FASTA.identifier(record))*" | SED = "* #sed should be changed in the future
                string(currminim)[1:5]*" | Pos = "*string(CMI+1)*":"*string(CMI+
                windowsize+1)*thrbuff,
                LongSequence{DNAAlphabet{4}}(seq[i-buff:i+windowsize-1+buff]))
                
                #return record to user
                if mode == "write"
                    FASTA.Writer(open(path, "a"), width = 95) do writer
                        write(writer, rec)
                    end
                elseif mode == "print"
                    println(">"*String(FASTX.FASTA.description(rec)))
                    println(getSeq(rec))
                else #if mode == "return"
                    push!(resultVec, rec)
                end

                #variable reset 
                currminim = currSqrEuc
                stop = true
            end
        end
    end
end

export gma

#looks like the main issue is that alot of the time is spent reading in the record...

"""
    eucGma(;k::Int64,
            record::FASTX.FASTA.Record,
            refVec::Vector{Float64},
            windowsize::Int64,
            kmerDict::Dict{LongSequence{DNAAlphabet{4}}, Int64},
            rv::Vector{Float64},
            resultVec::Vector{Float64} = Float64[],
            ScaleFactor::Float64 = 1.0) 

for devs only.

This version of the gma just puts the kmer dists into a vector instead of finding matches. 

I couldve just added an additional if statement in the actual gma but I wanted to optimize it as much as possible.
"""
function eucGma(; k::Int64,
                  record::FASTX.FASTA.Record,
                  refVec::Vector{Float64},
                  windowsize::Int64,
                  kmerDict::Dict{LongSequence{DNAAlphabet{4}}, Int64},
                  rv::Vector{Float64},
                  resultVec::Vector{Float64} = Float64[],
                  ScaleFactor::Float64 = 1.0) 

    sl = FASTX.FASTA.seqsize(record)
    if sl >= windowsize
        #initial operations for the first window 
        seq = getSeq(record) #getSeq is kinda slow. Ideally I wanna work with subseqs eventually. Even better is if I can read and edit at the same time
        curr = fasterKF(k,view(seq,1:windowsize),kmerDict,rv) #Theres an even faster way with unsigned ints and bits. I wonder if making so many temp arrays is a good idea
        currSqrEuc = Distances.sqeuclidean(refVec,curr)

        #initializing certain variables
        CMI = 2
        stop = true
        currminim = currSqrEuc

        for i in 1:(sl-windowsize) 
            #first & second operation
            #zeroth kmer
            @views zerokInt = kmerDict[view(seq,i:i+k-1)] #might be slightly faster to predefine
            @views currSqrEuc -= (refVec[zerokInt]-curr[zerokInt])^2 # a_old
            curr[zerokInt] -= 1
            @views currSqrEuc += (refVec[zerokInt]-curr[zerokInt])^2 # a_new

            #last kmer
            @views KendInt = kmerDict[view(seq,i+windowsize-k:i+windowsize-1)]
            @views currSqrEuc -= (refVec[KendInt]-curr[KendInt])^2 # b_old
            curr[KendInt] += 1
            @views currSqrEuc += (refVec[KendInt]-curr[KendInt])^2 # b_new

            #convert to kmer Distance and push to results
            push!(resultVec, currSqrEuc * ScaleFactor)
        end
    end
end

export eucGma

#using FlameGraphs, ProfileView, Profile
#Profile.clear(); @profile kmerGMA.writeQueryMatchN(6,LA,ref,d,190.0,289,50,"VicPacScan/vicpacscan.fasta"; rv = rV, thrbuff = "test")
#g = flamegraph()
#ProfileView.view()