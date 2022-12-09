#using FASTX, BioSequences, Distances, Main.KmerGMA, BenchmarkTools

#this entire script is a simpler version of the GMA that assumes 
const NUCLEOTIDE_BITS= Dict(DNA_A => unsigned(0),
                            DNA_C => unsigned(1),
                            DNA_G => unsigned(2),
                            DNA_T => unsigned(3))

const KmerType = Array{UInt32, 1}

#bins = zeros(eltype(KmerType), 4^k)
#mask = unsigned(4^k - 1)

#kmer counting in place, courtesy of the NextGenSeqUtils package
function kmer_count!(; str::LongSubSeq{DNAAlphabet{4}}, k::Int, 
    bins::KmerType, mask::UInt64, Nt_bits::Dict{DNA, UInt64})
    kmer = unsigned(0)
    for c in str[1:k-1]
        kmer = (kmer << 2) + Nt_bits[c]
    end
    for c in str[k:end]
        kmer = ((kmer << 2) & mask) + Nt_bits[c]
        bins[kmer + 1] += 1
    end
end

function kmer_count!(; str, k::Int, 
    bins::Vector, mask::UInt64, Nt_bits::Dict{DNA, UInt64})
    kmer = unsigned(0)
    for c in str[1:k-1]
        kmer = (kmer << 2) + Nt_bits[c]
    end
    for c in str[k:end]
        kmer = ((kmer << 2) & mask) + Nt_bits[c]
        bins[kmer + 1] += 1
    end
end

#fast reference generation
function gen_ref(path::String, k::Int, Nt_bits = NUCLEOTIDE_BITS)
    len = 0
    answer = zeros(Float64, 4^k)
    mask = unsigned(4^k - 1)
    open(FASTA.Reader,path) do reader
        for record in reader
            len += 1
            seq = getSeq(record)
            kmer_count!(str=seq,k=k,bins=answer,mask=mask,Nt_bits=Nt_bits)
        end
    end
    len = 1/len
    return answer.*len
end

#gen_ref("C:/Users/lu_41/.julia/dev/KmerGMA/test/Alp_V_ref.fasta",6)

"""
    gma_Nless(; k::Int64,
                record::FASTX.FASTA.Record,
                refVec::Vector{Float64},
                windowsize::Int64,
                thr::Float64,
                buff::Int64,
                curr_kmer_freq::Vector,
                thrbuff::String,
                mask::UInt64, 
                Nt_bits::Dict{DNA, UInt64},
                mode::String = "write",  #"return", "print"
                path::String = "noPath",
                resultVec::Vector{FASTA.Record} = FASTA.Record[],
                ScaleFactor::Float64 = 1.0)

Quick bit-based version of the gma (for devs only) that only works if no nucleotides are N

Does not rely on hashing and minimizes ooperations
"""
function gma_Nless(;
    k::Int64,
    record::FASTX.FASTA.Record,
    refVec::Vector{Float64},
    windowsize::Int64,
    thr::Float64,
    buff::Int64,
    curr_kmer_freq::Vector,
    thrbuff::String,
    mask::UInt64, 
    Nt_bits::Dict{DNA, UInt64},
    mode::String = "write",  #"return", "print"
    path::String = "noPath",
    resultVec::Vector{FASTA.Record} = FASTA.Record[],
    ScaleFactor::Float64 = 1.0)

    sequence_length = FASTX.FASTA.seqsize(record)
    if sequence_length >= windowsize
        #initial operations for the first window 
        seq = getSeq(record) #getSeq is kinda slow. Ideally I wanna work with subseqs eventually. Even better is if I can read and edit at the same time
        
        kmer_count!(str = view(seq,1:windowsize),k=k,
        bins=curr_kmer_freq,mask=mask,Nt_bits=Nt_bits)

        currSqrEuc = Distances.sqeuclidean(refVec,curr_kmer_freq)

        #initializing variables
        CMI = 2
        stop = true
        currminim = currSqrEuc

        #left kmer initialization
        left_kmer = unsigned(0)
        for c in seq[1:k-1]
            left_kmer = (left_kmer << 2) + Nt_bits[c]
        end

        #right kmer initialization
        right_kmer = unsigned(0)
        for c in seq[windowsize-k+1:windowsize-1]
            right_kmer = (right_kmer << 2) + Nt_bits[c]
        end

        for i in 1:(sequence_length-windowsize) #treat the first kmer as the -1th
            #first & second operation. I think a bloom filter might help with performance
            #zeroth kmer
            left_kmer = ((left_kmer << 2) & mask) + Nt_bits[seq[i+k-1]]
            left_ind = -~left_kmer
            @views currSqrEuc -= (refVec[left_ind]-curr_kmer_freq[left_ind])^2 # a_old
            curr_kmer_freq[left_ind] -= 1
            @views currSqrEuc += (refVec[left_ind]-curr_kmer_freq[left_ind])^2 # a_new

            #last kmer
            right_kmer = ((right_kmer << 2) & mask) + Nt_bits[seq[i+windowsize-1]]
            right_ind = -~right_kmer
            @views currSqrEuc -= (refVec[right_ind]-curr_kmer_freq[right_ind])^2 # b_old
            curr_kmer_freq[right_ind] += 1
            @views currSqrEuc += (refVec[right_ind]-curr_kmer_freq[right_ind])^2 # b_new

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
                #to get the record with the buffer, need to check if it exceeds the end of the sequence or is negative.
                left_buffer = i - buff
                if left_buffer < 0; left_buffer = 0 end 

                right_buffer = i + windowsize - 1 + buff 
                if right_buffer > sequence_length; right_buffer = sequence_length end  
                
                #create record 
                rec = FASTA.Record(String(FASTA.identifier(record))*" | SED = "* #sed should be changed in the future
                string(currminim)[1:5]*" | Pos = "*string(CMI+1)*":"*string(CMI+
                windowsize+1)*thrbuff,
                LongSequence{DNAAlphabet{4}}(seq[left_buffer:right_buffer]))
                
                #return record to user
                if mode == "write"
                    FASTA.Writer(open(path, "a"), width = 95) do writer # "a" is very important
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

export gma_Nless

"""
"""
function gma_Nless_API(; #FASTQ, RNA and AA compaibility will be added in the future. Also distance metric may be changed in future
    genome::Any, 
    ref::Any,
    thr::Float64, 
    Nt_bits = NUCLEOTIDE_BITS, #shouldnt even need this 
    mode::String = "return",
    fileloc::String = "noPath",
    k::Int64 = 6, 
    windowsize::Int64 = 0, 
    buffer::Int64 = 50,
    results::Vector{FASTA.Record} = FASTA.Record[])
 
    #managing undeclared variables
    k < 4 && println("Such a low k value may not yield the most accurate results")
    if windowsize == 0; windowsize = avgRecLen(ref) end
 
    #variables needed for the GMA
    RV = fill(0.0,4^k)
    ScaleFactor = 1/(2*k)
    threshold_buffer_tag = " | thr = "*string(round(thr))*" | buffer = "*string(buffer) # Maybe it would be wise to put which record it is
    mask = unsigned(4^k - 1)
 
    #handling KFD 
    refKFV = nothing
    if typeof(ref) == String
       refKFV = gen_ref(ref,k) #generation of kmer frequency dict
    else #assume its dict 
       refKFV = ref
    end
    
    #got rid of threshold finder for now. 20 is a good enough threshold 

    #genome mining, for string and vector of strings which is the database.
    
    if typeof(genome) == String 
        open(FASTA.Reader, genome) do io
            for rec in io
                fill!(RV,0.0)
                gma_Nless(k=k, record = rec, refVec = refKFV,
                windowsize = windowsize,
                path=fileloc, thr=thr, buff=buffer,
                curr_kmer_freq = RV, 
                thrbuff=threshold_buffer_tag,
                mask = mask, Nt_bits = Nt_bits,
                mode = mode, resultVec = results,
                ScaleFactor = ScaleFactor)
            end
        end
    else #for vector of strings
        for str in genome
            open(FASTA.Reader, str) do io
                for rec in io
                    fill!(RV,0.0)
                    gma_Nless(k=k, record = rec, refVec = refKFV,
                    windowsize = windowsize, 
                    path=fileloc, thr=thr, buff=buffer,
                    curr_kmer_freq = RV, 
                    thrbuff=threshold_buffer_tag,
                    mask = mask, Nt_bits = Nt_bits,
                    mode = mode, resultVec = results,
                    ScaleFactor = ScaleFactor)
                end
            end
        end
    end
    if mode == "return"; (return results) end
 end
 
 export GMA_Nless_API

 """
 refile = "C:/Users/lu_41/.julia/dev/KmerGMA/test/Alp_V_ref.fasta"
 genBank_long = "C:/Users/lu_41/.julia/dev/KmerGMA/test/Loci.fasta"
 genBank = "C:/Users/lu_41/.julia/dev/KmerGMA/test/Alp_V_locus.fasta"
 VicPac = "C:/Users/lu_41/Desktop/Sofo Prok/VicPac32.fna"

reference = gen_ref(refile,6)
@benchmark gma_Nless_API(genome=genBank,ref=reference,thr=30.0, windowsize=289)
#holy moly this speed is about 25 times faster 
"""
"""
BenchmarkTools.Trial: 4296 samples with 1 evaluation.
 Range (min … max):  642.300 μs …   5.791 ms  ┊ GC (min … max): 0.00% … 73.16%
 Time  (median):       1.230 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):     1.160 ms ± 362.735 μs  ┊ GC (mean ± σ):  1.13% ±  4.63%

                           ▆ █ ▆
  ▇▃█▃▇▃▃▂▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▄█▆█▆▅▄▃▃▃▃▃▂▂▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁ ▂
  642 μs           Histogram: frequency by time         1.94 ms <

 Memory estimate: 202.14 KiB, allocs estimate: 138.
"""

#@benchmark gma_Nless_API(genome=genBank_long,ref=reference,thr=30.0, windowsize=289)
"""
BenchmarkTools.Trial: 385 samples with 1 evaluation.
 Range (min … max):   7.033 ms … 17.993 ms  ┊ GC (min … max): 0.00% … 12.63%
 Time  (median):     13.912 ms              ┊ GC (median):    0.00%
 Time  (mean ± σ):   12.985 ms ±  2.702 ms  ┊ GC (mean ± σ):  0.53% ±  2.66%

                                          ▄█▆▂
  ▃▄▅▅▅▃▃▄▃▁▁▂▁▂▁▂▁▂▃▁▂▁▁▂▁▁▁▁▁▁▁▁▂▁▂▂▁▁▃▅████▆▅▅▄▃▄▃▄▃▃▃▃▃▂▂ ▃
  7.03 ms         Histogram: frequency by time        16.7 ms <

 Memory estimate: 1.43 MiB, allocs estimate: 335.
"""

#@benchmark findGenes()
