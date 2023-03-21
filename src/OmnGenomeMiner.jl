using BioAlignments, FASTX

# alot of the subprocesses can be written into functions for readability but for performance sake everything is kept in the function

function Omn_KmerGMA!(;
    genome_path::String,
    refVecs::Vector{Vector{Float64}},
    windowsizes::Vector{Int64}, 
    consensus_seqs::Vector{Seq},
    resultVec::Vector{FASTA.Record},
    k::Int64 = 6,
    ScaleFactor::Real = 0.16666666666666666666666666666666666666666666666667,
    mask::UInt64 = unsigned(4095), 
    thr_vec::Kfv = Float64[35,31,38,34,27,27], 
    buff::Int64 = 50, # if higher should increase accuracy
    Nt_bits::Dict{DNA, UInt64} = NUCLEOTIDE_BITS,
    score_model::AffineGapScoreModel{Int64} = AffineGapScoreModel(EDNAFULL, gap_open=-200, gap_extend=-1),
    genome_pos::Int64 = 0,
    align_hits::Bool = true,
    get_hit_loci::Bool = false,
    hit_loci_vec::Vector{Int} = Int[],
    get_aligns::Bool = false,
    align_vec::Vector{SubAlignResult} = SubAlignResult[],
    do_return_dists::Bool = false,
    dist_vec_vec::Vector{Vector{Float64}} = [Float64[] for _ in 1:6])

    # setup all variables
    len_KFVs::Int64 = length(windowsizes)
    kmerDist_vec = Float64[0 for _ in 1:len_KFVs]
    all_curr_kmer_freq = [zeros(4^k) for _ in 1:len_KFVs]
    curr_mins = Float64[Float64(100) for _ in 1:len_KFVs]
    CMIs = [1 for _ in 1:len_KFVs]
    stops = [true for _ in 1:len_KFVs]
    maxws = maximum(windowsizes)

    right_kmer_vec = UInt[unsigned(0) for _ in 1:len_KFVs]

    open(FASTX.FASTA.Reader, genome_path) do reader 
        for record in reader 
            # Operation on first window 
            seq::Seq = getSeq(record)
            sequence_length = length(seq)
            prev_hit_range = 0:0 # super important to avoid duplicate hits
            
            # initialize vars for the first window
            for ind in 1:len_KFVs
                if sequence_length < windowsizes[ind]; continue end #! need to sort out 
                
                fill!(all_curr_kmer_freq[ind], 0)

                kmer_count!(str = view(seq, 1:windowsizes[ind]), k = k,
                bins = all_curr_kmer_freq[ind], mask = mask)
                
                kmerDist_vec[ind] = curr_mins[ind] = (1/(2*k))*Distances.sqeuclidean(refVecs[ind], all_curr_kmer_freq[ind])

                CMIs[ind], stops[ind] = 1, true

                right_kmer_vec[ind] = unsigned(0)
                for c in view(seq, windowsizes[ind]-k+1:windowsizes[ind]-1) # here is the problem
                    right_kmer_vec[ind] = (right_kmer_vec[ind] << 2) + Nt_bits[c]
                end
            end

            left_kmer = unsigned(0)
            for c in view(seq, 1:k-1)
                left_kmer = (left_kmer << 2) + Nt_bits[c]
            end

            # iteration over the current record
            for i in 1:(sequence_length-maxws)

                # change left kmer
                left_kmer = ((left_kmer << 2) & mask) + Nt_bits[seq[i+k-1]]
                left_ind = 1+left_kmer

                for ind in 1:len_KFVs
                    # change right kmer
                    right_kmer_vec[ind] = ((right_kmer_vec[ind] << 2) & mask) + Nt_bits[seq[i+windowsizes[ind]-1]]
                    right_ind = 1 + right_kmer_vec[ind]

                    # simplified operation to update the distance
                    kmerDist_vec[ind] += ScaleFactor*(
                        all_curr_kmer_freq[ind][right_ind] + refVecs[ind][left_ind] -
                        refVecs[ind][right_ind] - all_curr_kmer_freq[ind][left_ind] + 1)

                    all_curr_kmer_freq[ind][left_ind] -= 1
                    all_curr_kmer_freq[ind][right_ind] += 1
                
                    kmerDist = kmerDist_vec[ind]
                    if do_return_dists; push!(dist_vec_vec[ind], kmerDist) end 

                    # minima finder (worse matches may sometimes be used but should be fine)
                    if kmerDist < thr_vec[ind]
                        if kmerDist < curr_mins[ind]
                            curr_mins[ind] = kmerDist
                            CMIs[ind] = i
                            stops[ind] = false
                        end
                    
                    # process actual hits
                    elseif !stops[ind]
                        stops[ind] = true
                        curr_mins[ind] = kmerDist 
                        CMI = CMIs[ind]

                        if !(CMI in prev_hit_range)
                            hit_left_ind, hit_right_ind = max(CMI-buff,1), min(CMI+windowsizes[ind]-1+buff,sequence_length)
                            seq_UnitRange = hit_left_ind:hit_right_ind

                            if align_hits
                                aligned_obj = pairalign(SemiGlobalAlignment(), consensus_seqs[ind], view(seq, seq_UnitRange), score_model)
                                if get_aligns; push!(align_vec, aligned_obj) end

                                aligned_UnitRange = cigar_to_UnitRange(aligned_obj)
                                seq_UnitRange = max(1, hit_left_ind+first(aligned_UnitRange)-1):min(hit_left_ind+last(aligned_UnitRange)-1, sequence_length)
                            end

                            # second overlap check
                            if last(seq_UnitRange) < first(prev_hit_range) || first(seq_UnitRange) > last(prev_hit_range)
                                if get_hit_loci; push!(hit_loci_vec, first(seq_UnitRange)+genome_pos) end

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
                                prev_hit_range = seq_UnitRange # !! #prev_hit_range = min(first(seq_UnitRange), first(prev_hit_range)):max(last(seq_UnitRange), last(prev_hit_range)) # could be even more/less conservative if wanted
                            end
                        end
                    end
                end
            end
            genome_pos += sequence_length
        end
    end
end

export Omn_KmerGMA!

#write_results(test_res, "new_testing.fasta") # took 6*50 = 300 seconds as expected, for 6 subsets, which is exactly 6 times longer than the O(n) version

# known edgecases not handled: 
# - the first window is a hit/close match
# - multiple hits under threshold - easy fix is to have new bound and reset minima finder, but might produce mroe false pos