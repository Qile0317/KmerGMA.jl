#using Pkg
#Pkg.add(PackageSpec(name="NextGenSeqUtils", rev="Missing-LongCharSeq-fix", url = "https://github.com/MurrellGroup/NextGenSeqUtils.jl.git"))

#this version just treates N as another nucleotide. it actually seems to work well as Ns incur a high sqeuclidean dist.
#EucVec version of GMA
"""
    eucGma(;
        k::Int64,
        record::FASTX.FASTA.Record,
        refVec::Vector{Float64},
        windowsize::Int64,
        kmerDict::Dict{LongSequence{DNAAlphabet{4}}, Int64},
        thr::Float64,
        buff::Int64,
        rv::Vector{Float64})

Testing version of the GMA that returns a euclidean vector instead of the matches.
"""
function eucGMA(;
    k::Int64,
    record::FASTX.FASTA.Record,
    refVec::Vector{Float64},
    windowsize::Int64,
    kmerDict::Dict{LongSequence{DNAAlphabet{4}}, Int64},
    thr::Float64,
    buff::Int64,
    rv::Vector{Float64})

    #initial operations
    seq = getSeq(record)
    @views curr = fasterKF(k,seq[1:windowsize],kmerDict,rv) #Theres an even faster way with unsigned ints and bits.
    currSqrEuc = Distances.sqeuclidean(refVec,curr)
    EucVec = [currSqrEuc]

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
        push!(EucVec, currSqrEuc)
    end
    return EucVec
end

export eucGMA

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

    #initial operations
    seq = getSeq(record)
    @views curr = fasterKF(k,seq[1:windowsize],kmerDict,rv) #Theres an even faster way with unsigned ints and bits.
    currSqrEuc = Distances.sqeuclidean(refVec,curr)

    #defining some variables needed later
    CMI = 2
    stop = true
    currminim = currSqrEuc

    for i in 1:(length(seq)-windowsize) 
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
        if currSqrEuc < thr
            if currSqrEuc < currminim
                currminim = currSqrEuc
                CMI = i
                stop = false
            end
        else
            if !stop #in future account for if buffer exceeds end or front
                #create the record of the match
                rec = FASTA.Record(String(FASTA.identifier(record))*" | SED = "*
                string(currminim)[1:5]*" | Pos = "*string(CMI+1)*":"*string(CMI+
                windowsize+1)*thrbuff,
                LongSequence{DNAAlphabet{4}}(seq[i-buff:i+windowsize-1+buff]))

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
function test_gma(;
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

    #initial operations
    seq = FASTA.sequence(LongSequence{DNAAlphabet{4}}, record)
    @views curr = fasterKF(k,seq[1:windowsize],kmerDict,rv) #Theres an even faster way with unsigned ints and bits.
    currSqrEuc = Distances.sqeuclidean(refVec,curr)

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
        if currSqrEuc < thr
            if currSqrEuc < currminim
                currminim = currSqrEuc
                CMI = i
                stop = false
            end
        else
            if stop == false #in future account for if buffer exceeds end or front
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
end
export test_gma

#using FlameGraphs, ProfileView, Profile
#Profile.clear(); @profile kmerGMA.writeQueryMatchN(6,LA,ref,d,190.0,289,50,"VicPacScan/vicpacscan.fasta"; rv = rV, thrbuff = "test")
#g = flamegraph()
#ProfileView.view()

"""
TO DO/TO ASK LIST:
    -> I gotta adjust for the new FASTX update bc substrs now are StringViews...
    -> BIG THING: In future I want to include D and J genes for comparison, and
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
