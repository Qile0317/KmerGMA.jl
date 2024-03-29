# unfinished version that is multithreadable
using BioSequences, FASTX, Distances, BioAlignments, Random, StaticArrays

# to account for multi-threading, a version for each record needs to be done and needs to use channels.
# Haven't quite yet figured out use of Channels, so this is just the naive solution to push to seperate vectors

# single_record gma for multithreading - unsure if inline works well w multithreading, need to benchmark and not sure if inlining is a good idea
@inline function record_KmerGMA!(;
    record::FASTX.FASTA.Record,
    refVec::SVector,
    curr_kmer_freq_vec::Vector{Vector{Int}}, # not sure if using Ref() could be better
    consensus_refseq::Seq,
    resultVec_vec::Vector{Vector{FASTA.Record}},
    k::Int64 = 6,
    windowsize::Int64 = 289,
    thr::Union{Int64, Float64} = 30,
    buff::Int64 = 50,
    mask::UInt64 = unsigned(4095), 
    Nt_bits::DnaBits = NUCLEOTIDE_BITS,
    ScaleFactor::Float64 = 0.166666666666666666666666666666666666666666666667,
    initial_scale_factor::Float64 = 0.0833333333333333333333333333333333333333333333,
    do_align::Bool = true,
    score_model::AffineGapScoreModel{Int64} = AffineGapScoreModel(EDNAFULL, gap_open=-69, gap_extend=-1))

    # get sequence
    seq::Seq = getSeq(record)
    sequence_length::Int = FASTX.FASTA.seqsize(record) # check if using 0.25*sizeof(seq) is better
    if sequence_length < windowsize; return end 
    thread_ind = Threads.threadid()

    #initial operations for the first window  
    @inbounds fill!(curr_kmer_freq_vec[thread_ind], 0) 
    @inbounds kmer_count!(str = view(seq, 1:windowsize),
        str_len = windowsize, k = k, bins = curr_kmer_freq_vec[thread_ind],
        mask = mask, Nt_bits = Nt_bits)

    @inbounds kmerDist::Float64 = initial_scale_factor *
        Distances.sqeuclidean(refVec, curr_kmer_freq_vec[thread_ind])

    left_kmer = unsigned(0); for i in 1:k-1
        @inbounds left_kmer = (left_kmer << 2) | Nt_bits[seq[i]]
    end

    right_kmer = unsigned(0); for i in windowsize-k+2:windowsize
        @inbounds right_kmer = (right_kmer << 2) | Nt_bits[seq[i]]
    end

    CMI::Int, stop::Bool, currminim::Float64, goal_ind::Int = 2, true, kmerDist, 0
    
    for (i_left, i_right) in zip(k:sequence_length-windowsize+k-1, windowsize+1:sequence_length)

        @inbounds left_kmer::UInt = ((left_kmer << 2) & mask) | Nt_bits[seq[i_left]]
        left_ind = left_kmer + 1

        @inbounds right_kmer::UInt = ((right_kmer << 2) & mask) | Nt_bits[seq[i_right]]
        right_ind = right_kmer + 1

        # simplified operation to update the distance
        if left_ind != right_ind
            @inbounds kmerDist += ScaleFactor * (1 +
                curr_kmer_freq_vec[thread_ind][right_ind] + refVec[left_ind] -
                refVec[right_ind] - curr_kmer_freq_vec[thread_ind][left_ind])
            
            # update the current kmer freq vector
            @inbounds curr_kmer_freq_vec[thread_ind][left_ind] -= 1
            @inbounds curr_kmer_freq_vec[thread_ind][right_ind] += 1
        end

        # minima finder
        if kmerDist < thr
            if kmerDist < currminim
                currminim = kmerDist
                CMI = i_left + 1 
                stop = false
            end
       
        # hit processing
        elseif !stop; stop = true
            if CMI > goal_ind; goal_ind = CMI + windowsize - 1

                seq_UnitRange = max(CMI - buff, 1):min(CMI + windowsize - 1 + buff, sequence_length)
                if do_align
                    seq_UnitRange = align_unitrange(
                        seq, seq_UnitRange, consensus_refseq, windowsize, sequence_length, score_model)
                end

                push!(resultVec_vec[Threads.threadid()], FASTA.Record(
                    FASTA.identifier(record)*
                        " | dist = "*string(round(currminim, digits = 2))*
                        " | MatchPos = $seq_UnitRange"*
                        " | Len = "*string(last(seq_UnitRange)-first(seq_UnitRange) + 1), 
                    view(seq, seq_UnitRange)
                ))
                currminim = kmerDist
            end
        end
    end
end

export record_KmerGMA!
"""
# function to take care of the actual threading from the API
function multi_threaded_KmerGMA!(;
    genome_path::String,
    refVec::Vector{Float64},
    consensus_refseq::Seq,
    resultVec_vec::Vector{Vector{FASTA.Record}},
    k::Int64 = 6,
    windowsize::Int64 = 289,
    thr::Union{Int64, Float64} = 30,
    buff::Int64 = 50,
    Nt_bits::DnaBits = NUCLEOTIDE_BITS,
    do_align::Bool = true,
    gap_open_score::Int = -69,
    gap_extend_score::Int = -1
    )

    num_kmers = 2 << ((2*k)-1)
    refVec = SVector{num_kmers}(refVec) # takes ab 6 ms
    curr_kmer_freq_vec = [zeros(Int, num_kmers) for _ in 1:Threads.nthreads()] # seems like MVector is surprisingly slower
    mask = unsigned(num_kmers - 1)
    ScaleFactor = 1/k
    initial_scale_factor = ScaleFactor * 0.5
    score_model = AffineGapScoreModel(EDNAFULL,
        gap_open = gap_open_score, gap_extend = gap_extend_score)

    open(FASTX.FASTA.Reader, genome_path) do reader 
        @sync for record in reader
            Threads.@spawn record_KmerGMA!(
                record = record, refVec = refVec,
                curr_kmer_freq_vec = curr_kmer_freq_vec,
                consensus_refseq = consensus_refseq, # does this create a copy? need to test if theres any improvements if i used subseq
                resultVec_vec = resultVec_vec, k = k,
                windowsize = windowsize, thr = thr, buff = buff,
                mask = mask, Nt_bits = Nt_bits, ScaleFactor = ScaleFactor,
                initial_scale_factor = initial_scale_factor,
                do_align = do_align, score_model = score_model
            )
        end
    end
end

#export threaded_KmerGMA!
"""
# so continually spawning threads destroys the memory, but a novel approach with channels is being implemented.