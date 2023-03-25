using BioSequences, FASTX, Distances, BioAlignments, Random, StaticArrays

# lots of memory is consumed because the records become stored in memory

@inline function ac_gma_testing!(;
    genome_path::String,
    refVec::Vector{Float64},
    consensus_refseq::Seq,
    k::Int64 = 6,
    windowsize::Int64 = 289,
    thr::Union{Int64, Float64} = 33.5,
    buff::Int64 = 50,
    mask::UInt64 = unsigned(4095), 
    Nt_bits::DnaBits = NUCLEOTIDE_BITS,
    ScaleFactor::Float64 = 0.1666666666666666666666666666666666666666666666667, # 1/k, NOT 1/2k !!
    do_align::Bool = true,
    score_model::AffineGapScoreModel{Int64} = AffineGapScoreModel(EDNAFULL, gap_open=-69, gap_extend=-1),
    do_return_dists::Bool = false, 
    do_return_align::Bool = false,
    get_hit_loci::Bool = false,
    dist_vec = Float64[], result_align_vec = [], hit_loci_vec = Int[],
    genome_pos::Int = 0, resultVec::Vector{FASTA.Record} = FASTA.Record[])

    # convert refVec to static vector for speed improvement
    refVec = SVector{2 << ((2*k)-1)}(refVec) # takes ab 6 secs
    curr_kmer_freq = zeros(Int, 2 << ((2*k)-1)) # seems like MVector is slower

    open(FASTX.FASTA.Reader, genome_path) do reader 
        for record in reader; seq::Seq = getSeq(record)

            sequence_length::Int = FASTX.FASTA.seqsize(record)
            if sequence_length < windowsize; continue end

            #initial operations for the first window  
            fill!(curr_kmer_freq, 0)
            kmer_count!(str = view(seq, 1:windowsize), k = k,
                bins = curr_kmer_freq, mask = mask, Nt_bits = Nt_bits)

            kmerDist::Float64 = ScaleFactor * 0.5 *
                Distances.sqeuclidean(refVec, curr_kmer_freq)

            left_kmer = unsigned(0); for nt in view(seq, 1:k-1)
                left_kmer = (left_kmer << 2) | Nt_bits[nt]
            end

            right_kmer = unsigned(0); for nt in view(seq, windowsize-k+1:windowsize-1)
                right_kmer = (right_kmer << 2) | Nt_bits[nt]
            end

            CMI, stop, currminim, i, goal_ind = 2, true, kmerDist, 0, 0
            
            for nt in view(seq, k:sequence_length-windowsize+1); i += 1
                
                left_kmer::UInt = ((left_kmer << 2) & mask) | Nt_bits[nt]
                left_ind = left_kmer + 1

                @inbounds right_kmer::UInt = ((right_kmer << 2) & mask) | Nt_bits[seq[i+windowsize-1]]
                right_ind = right_kmer + 1

                # simplified operation to update the distance - it is right, right??
                @inbounds kmerDist += ScaleFactor * (1 +
                    curr_kmer_freq[right_ind] + refVec[left_ind] -
                    refVec[right_ind] - curr_kmer_freq[left_ind])

                if do_return_dists; push!(dist_vec, kmerDist) end 

                # update the current kmer freq vector
                @inbounds curr_kmer_freq[left_ind] -= 1
                @inbounds curr_kmer_freq[right_ind] += 1

                # minima finder
                if kmerDist < thr
                    if kmerDist < currminim
                        currminim = kmerDist
                        CMI = i
                        stop = false
                    end
                
                # hit processing 
                elseif !stop
                    if CMI > goal_ind; goal_ind = CMI + windowsize - 1
                        # alignment and matched unitrange
                        left_ind, right_ind = max(CMI-buff,1), min(CMI+windowsize-1+buff,sequence_length)

                        if do_align 
                            aligned_obj = pairalign(SemiGlobalAlignment(),consensus_refseq,view(seq,left_ind:right_ind),score_model)
                            if do_return_align; push!(result_align_vec, aligned_obj) end
                            aligned_UnitRange = cigar_to_UnitRange(aligned_obj)
                            seq_UnitRange = max(1, left_ind+first(aligned_UnitRange)-1):min(left_ind+last(aligned_UnitRange)-1, sequence_length)
                        else
                            seq_UnitRange = left_ind:right_ind
                        end

                        #create record
                        rec = FASTA.Record(
                            FASTA.identifier(record)*
                                " | dist = "*string(round(currminim, digits = 2))*
                                " | MatchPos = $seq_UnitRange"*
                                " | GenomePos = $genome_pos"*
                                " | Len = "*string(last(seq_UnitRange)-first(seq_UnitRange)),
                            view(seq, seq_UnitRange)
                        )

                        if get_hit_loci; push!(hit_loci_vec, first(seq_UnitRange)+genome_pos) end 
                        push!(resultVec, rec)
                        currminim = kmerDist
                        stop = true
                    end
                end
            end
            genome_pos += sequence_length
        end
    end
end

export ac_gma_testing!

# should use profiler for optimization