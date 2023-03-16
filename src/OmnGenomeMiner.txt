# alternative approach for multiple RVs that account for duplicates. Unfinished and Bugged at the moment.

function Omn_KmerGMA!(;
    genome_path::String,
    refVecs::Vector{Vector{Float64}},
    k::Int64,
    ScaleFactor::Real,
    windowsizes::Vector{Int64}, 
    mask::UInt64, 
    consensus_seqs::Vector{Seq},
    thr::Union{Int64, Float64} = 30.0, # could be a vec
    buff::Int64 = 50, # could be a vec 
    Nt_bits::Dict{DNA, UInt64} = NUCLEOTIDE_BITS,
    score_model::AffineGapScoreModel{Int64} = AffineGapScoreModel(EDNAFULL, gap_open=-69, gap_extend=-1),
    align_hits::Bool = true,
    do_return_dists::Bool = false)

    curr_kmer_freq::Vector{Vector{Float64}} = [zeros(4096) for _ in 1:3]
    resultVec = FASTA.Record[]
    len_KFVs::Int64 = length(refVecs)
    currSqrEuc_vec = Float64[0 for _ in 1:len_KFVs]
    genome_pos::Int64 = 0
    if align_hits; result_align_vec = [] end
    if do_return_dists; dist_vec_vec = [Float64[] for _ in 1:len_KFVs] end 

    CMIs = [1 for _ in 1:len_KFVs]
    stops = [true for _ in 1:len_KFVs]
    maxws = maximum(windowsizes)

    # begin iteration (could be made into a function)
    open(FASTX.FASTA.Reader, genome_path) do reader 
        for record in reader 
            # Operation on first window 
            seq::LongSequence{DNAAlphabet{4}} = KmerGMA.getSeq(record)
            sequence_length = length(seq)
            proceed::Bool, goal_ind = true, 0

            for ind in 1:len_KFVs
                if sequence_length < windowsizes[ind]; continue end #! need to sort out 

                fill!(curr_kmer_freq[ind], 0)

                KmerGMA.kmer_count!(str = view(seq, 1:windowsizes[ind]), k = k,
                bins = curr_kmer_freq[ind], mask = mask, Nt_bits = Nt_bits)
                
                currSqrEuc_vec[ind] =  Distances.sqeuclidean(refVecs[ind], curr_kmer_freq[ind])
            end

            #initializing variables for distance
            fill!(CMIs, 1); fill!(stops, true)
            curr_mins = [num for num in currSqrEuc_vec]

            left_kmer = unsigned(0)
            for c in view(seq, 1:k-1)
                left_kmer = (left_kmer << 2) + Nt_bits[c]
            end

            right_kmer_vec = UInt[unsigned(0) for _ in 1:len_KFVs]
            for i in 1:len_KFVs
                for c in view(seq, windowsizes[i]-k+1:windowsizes[i]-1) # here is the problem
                    right_kmer_vec[i] = (right_kmer_vec[i] << 2) + Nt_bits[c]
                end
            end
            right_ind_vec = UInt[unsigned(0) for _ in 1:len_KFVs]

            # iteration 
            for i in 1:(sequence_length-maxws) # rip, assumes even distribution

                left_kmer = ((left_kmer << 2) & mask) + Nt_bits[seq[i+k-1]]
                left_ind = left_kmer + 1 
                #if left_ind > 4000; print(" | $left_ind "); print(seq[i+k-1]) end end

                for ind in 1:len_KFVs
                    # change right kmer
                    right_kmer_vec[ind] = ((right_kmer_vec[ind] << 2) & mask) + Nt_bits[seq[i+windowsizes[ind]-1]]
                    right_ind_vec[ind] = right_kmer_vec[ind] + 1

                    # update left kmer
                    currSqrEuc_vec[ind] -= (refVecs[ind][left_ind] - curr_kmer_freq[ind][left_ind])^2 # problem, the left_ind is not right
                    currSqrEuc_vec[ind] += (refVecs[ind][left_ind] - curr_kmer_freq[ind][left_ind] + 1)^2

                    # update right kmer
                    currSqrEuc_vec[ind] -= (refVecs[ind][right_ind_vec[ind]] - curr_kmer_freq[ind][right_ind_vec[ind]])^2
                    currSqrEuc_vec[ind] += (refVecs[ind][right_ind_vec[ind]] - curr_kmer_freq[ind][right_ind_vec[ind]] - 1)^2
                    
                    curr_kmer_freq[ind][left_ind] -= 1
                    curr_kmer_freq[ind][right_ind_vec[ind]] += 1
                
                    kmerDist = currSqrEuc_vec[ind] * ScaleFactor 
                    if do_return_dists; push!(dist_vec_vec[ind], kmerDist) end

                    # minima finder (worse matches may sometimes be used, needs fix)
                    if kmerDist < thr
                        if kmerDist < curr_mins[ind]
                            curr_mins[ind] = kmerDist
                            CMIs[ind] = i
                            stops[ind] = false
                        end

                    elseif !stops[ind]
                        CMI = CMIs[ind]

                        if !proceed 
                            if CMI > goal_ind # problem; edge case is when theres another close... but ehh...
                                proceed = true #universal!
                            end
                        else 
                            goal_ind = CMI + windowsizes[ind] - 1
                            proceed = false 

                            # alignment and matched unitrange
                            left_ind, right_ind = max(CMI-buff,1), min(CMI+windowsizes[ind]-1+buff,sequence_length)

                            if align_hits
                                aligned_obj = pairalign(SemiGlobalAlignment(), consensus_seqs[ind], view(seq, left_ind:right_ind), score_model)
                                push!(result_align_vec, aligned_obj)
                                aligned_UnitRange = cigar_to_UnitRange(cigar(aligned_obj.aln.a.aln))
                                seq_UnitRange = max(1, left_ind+first(aligned_UnitRange)-1):min(left_ind+last(aligned_UnitRange)-1, sequence_length)
                            else
                                seq_UnitRange = left_ind:right_ind
                            end

                            #create and push record
                            push!(resultVec, FASTA.Record(
                                FASTA.identifier(record)*
                                    " | Dist = "*string(round(curr_mins[ind], digits = 2))*
                                    " | KFV = $ind"*
                                    " | MatchPos = $seq_UnitRange"*
                                    " | GenomePos = $genome_pos"*
                                    " | Len = "*string(last(seq_UnitRange)-first(seq_UnitRange)),
                                view(seq, seq_UnitRange)
                            ))

                            curr_mins[ind] = currSqrEuc_vec[ind]
                            stops .= true # dont think will help but :/
                        end
                    end
                end
            end
            genome_pos += sequence_length
        end
    if align_hits
        if do_return_dists; return resultVec, result_align_vec, dist_vec_vec end
        return resultVec, result_align_vec
    end
    if do_return_dists
        return resultVec, dist_vec_vec
    end
    return resultVec
end