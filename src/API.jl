"""
   findGenes(; genome::FASTX.FASTA.Reader,
               ref::FASTX.FASTA.Reader,
               fileloc::String,
               k::Int64 = 6,
               windowsize::Int64 = 0,
               print::Bool = false,
               thr::Int64 = 0,
               buffer::Int64 = 50,
               BLAST = true)

The main API to find all matches. Its impportant that the genome to be scanned is a reader object.

The ref argument is a reader of reference FASTA sequences. fileloc is simply the location of a blank file to read the matches to.

Unfinished. The print and BLAST arguments are currently useless.
"""
function findGenes(;genome::FASTX.FASTA.Reader,
   ref::FASTX.FASTA.Reader, fileloc::String,
   k::Int64 = 6, windowsize::Int64 = 0,
   print::Bool = false, thr::Int64 = 0,
   buffer::Int64 = 50, BLAST = true)

   #managing undeclared variables
   k < 4 && error("try a higher value like 6. It is most likely more accurate") #adressing k-value
   windowsize == 0 && (windowsize = avgRecLen(ref,true)) #finding of adequate windowsize

   #variables needed for the GMA
   KD = genKmers(k;withN=true) #generation of kmer dictionary
   refKFV = kfv(genRef(k,ref,KD),KD)
   RV = fill(0.0,5^k)
   threshold_buffer_tag = " | thr = "*string(thr)*" | buffer = "*string(buffer)
   thr == 0 && (thr = Float64(findthr(ref,refKFV,KD))) #SED threshold estimation

    #genome mining
    for record in genome
       gma(k=k, record = record, refVec =refKFV,
       windowsize = windowsize, kmerDict = KD,
       path=fileloc, thr=thr, buff=buffer,
       rv=copy(RV), #from benchmarking it seems copying isnt that slow?
       thrbuff=threshold_buffer_tag)
    end
    close(ref)
    close(genome)
    #if BLAST
      #blastseq = eqBLAST(output) #blasting
      #return blastseq #unfinished
   #end
end

export findGenes

#testing version, same code but doesnt write to a file and instead pushes to a vector
function test_gma(; k::Int64,
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
    @views curr = fasterKF(k,seq[1:windowsize],kmerDict,rv) #Theres an even faster way with unsigned ints and bits
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
                rec = FASTA.Record(FASTA.identifier(record),
                "| SED = "*string(currminim)[1:5]*" | Pos = "*string(CMI+1)*":"*string(CMI+windowsize+1)*thrbuff,
                seq[i-buff:i+windowsize-1+buff])

                #write in the record to the file
                push!(path, rec)

                #reset
                currminim = currSqrEuc
                stop = true
            end
        end
    end
end

open(FASTX.FASTA.Reader, "test/Alp_V_ref.fasta") do reference
    open(FASTX.FASTA.Reader, "test/Alp_V_locus.fasta") do target
        KD = genKmers(6,withN=true)
        k = 6
        refKFD = genRef(k,reference,KD) #generation of kmer frequency dict
        refKFV = kfv(refKFD,KD)
        inp = FASTA.Record[]

        test_gma(k=k,record=first(target),
        refVec = refKFV, windowsize = 289,
        kmerDict =KD,
        path = inp,
        thr = 250.0, buff = 20,
        rv= fill(0.0,5^6),
        thrbuff = " test ")

        print(inp)
    end
end


function testFindGenes(; genome::FASTX.FASTA.Reader,
    ref::FASTX.FASTA.Reader,
    k::Int64 = 6, windowsize::Int64 = 0,
    print::Bool = false, thr::Int64 = 0,
    buffer::Int64 = 50, BLAST = true)

    k < 4 && error("try a higher value like 6. It is most likely more accurate") #adressing k-value
    windowsize == 0 && (windowsize = avgRecLen(ref,true)) #finding of adequate windowsize

    #variables needed for the GMA
    KD = genKmers(k;withN=true) #generation of kmer dictionary
    refKFD = genRef(k,ref,KD) #generation of kmer frequency dict
    refKFV = kfv(refKFD,KD) #generation of kmer frequency dict
    RV = fill(0.0,5^k)
    threshold_buffer_tag = " | thr = "*string(thr)*" | buffer = "*string(buffer)
    results = FASTX.FASTA.Record[]
    thr == 0 && (thr = Float64(findthr(ref,refKFD,KD))) #SED threshold estimation

    #genome mining, code copied from gma()
    for rec in genome
        test_gma(k=k, record = rec, refVec =refKFV,
        windowsize = windowsize, kmerDict = KD,
        path=results, thr=thr, buff=buffer,
        rv=copy(RV), #from benchmarking it seems copying isnt that slow?
        thrbuff=threshold_buffer_tag)
    end
    close(ref)
    close(genome)
    return results
end
