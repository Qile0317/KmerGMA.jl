# unfortunately it seems like useing strobes dont really offer any advantage to kmers
# feels more promising if vectors were weighed following an analysis of the reference sequences

function StrobeGMA!(;
    genome_path::String,
    refVec::SVector, 
    consensus_refseq::Seq,

    s::Int = 2, w_min::Int = 3, w_max::Int = 5, q::Int = 5,
    windowsize::Int64 = 289,
    thr::Union{Int64, Float64} = 33.5,
    ScaleFactor::Float64 = 0.166666666666666666666666666667,

    buff::Int64 = 50,
    do_align::Bool = true,
    score_model::AffineGapScoreModel{Int64} = AffineGapScoreModel(EDNAFULL, gap_open=-69, gap_extend=-1),
    do_return_dists::Bool = false, # dangerous for memory if true
    do_return_align::Bool = false,
    get_hit_loci::Bool = false,

    dist_vec = Float64[], result_align_vec = [], hit_loci_vec = Int[],
    genome_pos::Int = 0, resultVec::Vector{FASTA.Record} = FASTA.Record[])

    curr_strobemer_freq::Vector{Float64} = zeros(2 << ((4*s)-1))
    k = w_max+s-1
    open(FASTX.FASTA.Reader, genome_path) do reader 
        for record in reader 

            seq::Seq = getSeq(record)
            goal_ind::Int64 = 0

            # length and edge case
            sequence_length = FASTX.FASTA.seqsize(record)
            if sequence_length < windowsize; continue end

            #initial operations for the first window  
            fill!(curr_strobemer_freq, 0)
            ungapped_strobe_2_mer_count!(view(seq, 1:windowsize), curr_strobemer_freq,
                s, w_min, w_max, q)
            #display(plot(curr_strobemer_freq))

            kmerDist = (1/(2*k))*Distances.sqeuclidean(refVec, curr_strobemer_freq)

            #initializing variables
            CMI, stop, currminim = 2, true, kmerDist

            for i in 1:(sequence_length-windowsize-1)
                # i-1 th kmer 
                @inbounds left_strobemer = get_strobe_2_mer(view(seq, i:i+k-1),
                    s, w_min, w_max, q; withGap = false)
                left_ind = as_UInt(left_strobemer) + 1
                
                # windowsize + 1th kmer 
                @inbounds right_strobemer = get_strobe_2_mer(view(seq, i+windowsize-k:i+windowsize),
                    s, w_min, w_max, q; withGap = false)
                right_ind = as_UInt(right_strobemer) + 1

                if left_ind != right_ind
                    @inbounds kmerDist += ScaleFactor*(1 +  
                        curr_strobemer_freq[right_ind] + refVec[left_ind] -
                        refVec[right_ind] - curr_strobemer_freq[left_ind])

                    # update the current kmer freq vector
                    @inbounds curr_strobemer_freq[left_ind] -= 1
                    @inbounds curr_strobemer_freq[right_ind] += 1
                end

                if do_return_dists; push!(dist_vec, kmerDist) end 

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

"""
    Strobemer_findGenes(;
        genome_path::String,
        ref_path::String,
        s::Int = 2,
        w_min::Int = 3,
        w_max::Int = 5,
        q::Int = 5,
        KmerDistThr::Union{Int64, Float64} = 30,
        buffer::Int64 = 50,
        do_align::Bool = true,
        do_return_dists::Bool = false,
        do_return_hit_loci::Bool = false,
        do_return_align::Bool = false,
        verbose::Bool = true)

An experimental homology searcher that uses Strobemers, specifically randstrobes with 2 sub-kmers. 
Very unoptimized and there is not distance threshold estimation yet. More documentation is to come.
"""
function Strobemer_findGenes(; genome_path::String, ref_path::String,
    s::Int = 2, w_min::Int = 3, w_max::Int = 5, q::Int = 5,
    KmerDistThr::Union{Int64, Float64} = 30, buffer::Int64 = 50, do_align::Bool = true,
    do_return_dists::Bool = false, do_return_hit_loci::Bool = false, do_return_align::Bool = false,
    verbose::Bool = true)

    RV, windowsize, consensus_refseq = gen_ref_ws_cons(ref_path; 
        s=s, w_min = w_min, w_max = w_max, q = q)
    
    hit_vector = FASTX.FASTA.Record[]
    dist_vec, hit_loci_vec, alignment_vec = Float64[], Int[], []  
    cumulative_length_in_genome = 0

    if verbose; @info "initializing iteration..." end 
    StrobeGMA!(
        genome_path = genome_path, refVec = RV, consensus_refseq = consensus_refseq,
        s = s, w_min = w_min, w_max = w_max, q = q,
        windowsize = windowsize, thr = KmerDistThr, ScaleFactor = 1/(w_max+s-1),
        buff = buffer,

        do_align = do_align, do_return_dists = do_return_dists,
        do_return_align = do_return_align, get_hit_loci = do_return_hit_loci,

        dist_vec = dist_vec,
        result_align_vec = alignment_vec,
        hit_loci_vec = hit_loci_vec, genome_pos = cumulative_length_in_genome,
        resultVec = hit_vector)

    if verbose; info_str = "genome mining completed successfully, returning vector of: vector of hits" end 
    output_vector = Any[hit_vector]
    if do_return_hit_loci; push!(output_vector, hit_loci_vec); if verbose; info_str *= ", vector of hit locations" end end
    if do_return_align; push!(output_vector, alignment_vec); if verbose; info_str *= ", vector of alignments" end end
    if do_return_dists; push!(output_vector, dist_vec) ; if verbose; info_str *= ", vector of kmer distances along the genome"end end 
    
    if verbose; @info info_str end
    return output_vector
end

export Strobemer_findGenes

"""
a = sgma(genome_path = "test/Loci.fasta", ref_path = "test/Alp_V_ref.fasta",
    s = 2, w_min = 3, w_max = 5, q = 5, KmerDistThr = 16,
    do_return_dists = true)

plot!(a[2])

res= FASTA.Record[]
distvec = Float64[]
refVec, ws, cons = gen_ref_ws_cons("test/Alp_V_ref.fasta", 6)
@time ac_gma_testing!(genome_path = "test/Loci.fasta",
    refVec = refVec, consensus_refseq = cons, resultVec = res,
    do_return_dists = true, dist_vec = distvec, do_align = false)

plot!(distvec)
"""