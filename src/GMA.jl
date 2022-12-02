#using Pkg
#Pkg.add(PackageSpec(name="NextGenSeqUtils", rev="Missing-LongCharSeq-fix", url = "https://github.com/MurrellGroup/NextGenSeqUtils.jl.git"))

#this version just treates N as another nucleotide. it actually seems to work well as Ns incur a high sqeuclidean dist.

####################################################################
#faster one for record, essential for GMA!!!!
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
          resultVec::Vector{FASTA.Record} = FASTA.Record[])

the main genome mining algorithm for a single record. either prints, or edit a vector or file in place

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
    resultVec::Vector{FASTA.Record} = FASTA.Record[])

    sl = FASTX.FASTA.seqsize(record)
    if sl >= windowsize
        #initial operations for the first window 
        seq = getSeq(record) #getSeq is kinda slow. Ideally I wanna work with subseqs eventually. Even better is if I can read and edit at the same time
        curr = fasterKF(k,view(seq,1:windowsize),kmerDict,rv) #Theres an even faster way with unsigned ints and bits.
        currSqrEuc = Distances.sqeuclidean(refVec,curr)

        #defining some variables needed later
        CMI = 2
        stop = true
        currminim = currSqrEuc

        for i in 1:(sl-windowsize) 
            #first & second operation
            #zeroth kmer
            @views zerokInt = kmerDict[view(seq,i:i+k-1)] #might be slightly faster to predefine
            @views currSqrEuc -= (refVec[zerokInt]-curr[zerokInt])^2 #MinusOld
            curr[zerokInt] -= 1
            @views currSqrEuc += (refVec[zerokInt]-curr[zerokInt])^2 #MinusNew

            #last kmer
            @views KendInt = kmerDict[view(seq,i+windowsize-k:i+windowsize-1)]
            @views currSqrEuc -= (refVec[KendInt]-curr[KendInt])^2 #plusOld
            curr[KendInt] += 1
            @views currSqrEuc += (refVec[KendInt]-curr[KendInt])^2 #PlusNew

            #third operation - the block below !stop's speed is irrelevant as its rare
            if currSqrEuc < thr
                if currSqrEuc < currminim
                    currminim = currSqrEuc
                    CMI = i
                    stop = false
                end
            elseif !stop #in future account for if buffer exceeds end or front
                #create the record of the match
                rec = FASTA.Record(String(FASTA.identifier(record))*" | SED = "*
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

                #reset
                currminim = currSqrEuc
                stop = true
            end
        end
    end
end

export gma

#looks like the main issue is that alot of the time is spent reading in the record...
#now that I think about it... what if I just use the manhattan distance?? wont it just be faster

#testing version that pushes instead of writing. Also it could actually be a viable alternative

#testing version, same code as gma but doesnt write to a file and instead pushes to a vector
#it only pushes nucleotides 10 to 20 for testing purposes
"""
    test_gma(;
        k::Int64,
        record::FASTX.FASTA.Record,
        refVec::Vector{Float64},
        windowsize::Int64,
        kmerDict::Dict{LongSequence{DNAAlphabet{4}}, Int64},
        path::Vector{FASTX.FASTA.Record},
        thr::Float64,
        buff::Int64,
        rv::Vector{Float64},
        thrbuff::String)

Testing version of the gma, it only returns [10:20] nucleotides of the matches
"""
function test_gma(; #this version is actually slightly worse than the GMA now 
    k::Int64,
    record::FASTX.FASTA.Record,
    refVec::Vector{Float64},
    windowsize::Int64,
    kmerDict::Dict{LongSequence{DNAAlphabet{4}}, Int64},
    thr::Float64,
    buff::Int64,
    rv::Vector{Float64},
    thrbuff::String = "placeholder",
    path::Vector{FASTX.FASTA.Record} = FASTA.Record[],
    euc::Bool = false)

    #initial operations
    seq = FASTA.sequence(LongSequence{DNAAlphabet{4}}, record)
    curr = fasterKF(k,view(seq, 1:windowsize),kmerDict,rv) #Theres an even faster way with unsigned ints and bits.
    currSqrEuc = Distances.sqeuclidean(refVec,curr)
    eucVec = [currSqrEuc]

    #defining some variables needed later
    CMI = 2
    stop = true
    currminim = currSqrEuc

    for i in 1:(length(seq)-windowsize) #big change: starts from 1 so i dont need to -1 nymore
        #first operation
        #zeroth kmer
        @views zerokInt = kmerDict[seq[i:i+k-1]] #might be slightly faster to predefine
        @views MinusOld = (refVec[zerokInt]-curr[zerokInt])^2
        curr[zerokInt] -= 1
        @views MinusNew = (refVec[zerokInt]-curr[zerokInt])^2

        #last kmer
        @views KendInt = kmerDict[seq[i+windowsize-k:i+windowsize-1]]
        @views PlusOld = (refVec[KendInt]-curr[KendInt])^2
        curr[KendInt] += 1
        @views PlusNew = (refVec[KendInt]-curr[KendInt])^2

        #second operation.
        currSqrEuc += PlusNew + MinusNew - PlusOld - MinusOld

        #third operation 
        if euc
            push!(eucVec, currSqrEuc)
        else 
            if currSqrEuc < thr
                if currSqrEuc < currminim
                    currminim = currSqrEuc
                    CMI = i
                    stop = false
                end
            elseif !stop
                #in future account for if buffer exceeds end or front
                #create the record of the match
                rec = FASTA.Record(String(FASTA.identifier(record))*" | SED = "*
                string(currminim)[1:5]*" | Pos = "*string(CMI+1)*":"*string(CMI+
                windowsize+1)*thrbuff,
                LongSequence{DNAAlphabet{4}}(seq[i-buff:i+windowsize-1+buff])[10:20])

                #write in the record to the file
                push!(path, rec)

                #reset
                currminim = currSqrEuc
                stop = true
            end
        end
    end
    if euc; (return eucVec) end
end
export test_gma

#using FlameGraphs, ProfileView, Profile
#Profile.clear(); @profile kmerGMA.writeQueryMatchN(6,LA,ref,d,190.0,289,50,"VicPacScan/vicpacscan.fasta"; rv = rV, thrbuff = "test")
#g = flamegraph()
#ProfileView.view()