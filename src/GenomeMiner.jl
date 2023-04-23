using BioSequences, FASTX, Distances, BioAlignments, Random, StaticArrays

# alot of the code would be cleaner if classes were used but classes are unfortunately very slow
function ac_gma_testing!(;
    genome_path::String,
    refVec::Vector{Float64},
    consensus_refseq::Seq,
    k::Int64 = 6,
    windowsize::Int64 = 289,
    thr::Union{Int64, Float64} = 33.5,
    buff::Int64 = 50,
    mask::UInt64 = unsigned(4095), 
    Nt_bits::DnaBits = NUCLEOTIDE_BITS,
    ScaleFactor::Float64 = 0.166666666666666666666666666666666666666666666667,
    do_align::Bool = true,
    result_align_vec::Vector{AlignResult} = AlignResult[],
    gap_open_score::Int = -69,
    gap_extend_score::Int = -1,
    do_return_dists::Bool = false, dist_vec = Float64[],
    do_return_align::Bool = false,
    get_hit_loci::Bool = false,
    hit_loci_vec::Vector{Int} = Int[],
    resultVec::Vector{FASTA.Record} = FASTA.Record[])
    
    genome_pos::Int = 0
    refVec = SVector{2 << ((2*k)-1)}(refVec) 
    curr_kmer_freq = zeros(Int, 2 << ((2*k)-1)) # seems like MVector is slower?
    score_model = AffineGapScoreModel(EDNAFULL, gap_open = gap_open_score, gap_extend = gap_extend_score)
    initial_scale_factor::Float64 = ScaleFactor * 0.5
    
    open(FASTX.FASTA.Reader, genome_path) do reader 
        for record in reader

            sequence_length = FASTX.FASTA.seqsize(record)
            seq::Seq = getSeq(record)
            
            if sequence_length < windowsize
                continue
            end 
            
            #initial operations for the first window  
            fill!(curr_kmer_freq, 0)
            kmer_count!(str = view(seq, 1:windowsize), str_len = windowsize,
                k = k, bins = curr_kmer_freq, mask = mask, Nt_bits = Nt_bits)

            kmerDist::Float64 = initial_scale_factor *
                Distances.sqeuclidean(refVec, curr_kmer_freq)

            left_kmer = unsigned(0); for i in 1:k-1
                @inbounds left_kmer = (left_kmer << 2) | Nt_bits[seq[i]]
            end

            right_kmer = unsigned(0); for i in windowsize-k+2:windowsize
                @inbounds right_kmer = (right_kmer << 2) | Nt_bits[seq[i]]
            end

            CMI::Int, stop::Bool, currminim::Float64, goal_ind::Int = 2, true, kmerDist, 0
            
            # iterative phase
            for (i_left, i_right) in zip(k:sequence_length-windowsize+k-1, windowsize+1:sequence_length)

                @inbounds left_kmer::UInt = ((left_kmer << 2) & mask) | Nt_bits[seq[i_left]]
                left_ind = left_kmer + 1

                @inbounds right_kmer::UInt = ((right_kmer << 2) & mask) | Nt_bits[seq[i_right]]
                right_ind = right_kmer + 1

                # simplified operation to update the distance
                if left_ind != right_ind
                    @inbounds kmerDist += ScaleFactor * (1 +
                        curr_kmer_freq[right_ind] + refVec[left_ind] -
                        refVec[right_ind] - curr_kmer_freq[left_ind])
                    
                    # update the current kmer freq vector
                    @inbounds curr_kmer_freq[left_ind] -= 1
                    @inbounds curr_kmer_freq[right_ind] += 1
                end

                if do_return_dists; push!(dist_vec, kmerDist) end 

                # minima finder
                if kmerDist < thr
                    if kmerDist < currminim
                        currminim = kmerDist
                        CMI = i_left
                        stop = false
                    end
               
                # hit processing 
                elseif !stop
                    stop = true
                    CMI += 1
                    if CMI > goal_ind
                        goal_ind = CMI + windowsize - 1
                        seq_UnitRange = (max(CMI - buff, 1)):(min(CMI + windowsize - 1 + buff, sequence_length))
                        if do_align
                            seq_UnitRange = align_unitrange(
                                seq, seq_UnitRange, consensus_refseq, windowsize, sequence_length, score_model, do_return_align, result_align_vec)
                        end
                        append_hit!(resultVec, record, seq, false, 0, currminim, seq_UnitRange, genome_pos)
                        if get_hit_loci; push!(hit_loci_vec, first(seq_UnitRange) + genome_pos) end
                        currminim = kmerDist
                    end
                end
            end
            genome_pos += sequence_length
        end
    end
end

export ac_gma_testing!

# should use profiler for optimization
"""
using BenchmarkTools
RV, ws, cons_seq = gen_ref_ws_cons(tf, 6)
res = FASTA.Record[]
@benchmark ac_gma_testing!(genome_path = test_mini_genome,
    refVec = RV, consensus_refseq = cons_seq,
    windowsize = ws, thr = 30, do_align = false, #do_overlap = false,
    resultVec = res)
    """ 
# looks like appending sequences for the overlap was super slow as expected