"""
# unfinished script to multithread 
using BioSequences, FASTX, Distances, BioAlignments, Random, StaticArrays

# lots of memory is consumed because the records become stored in memory
# gotta check the indexing for the left and right kmers...

# to account for multi-threading, a version for each record needs to be done and needs to use channels.
# is it a good idea to keep the original version 

# single_record gma for multithreading - unsure if inline works well w multithreading, need to benchmark
@inline function record_KmerGMA!(;
    record::FASTX.FASTA.Record,
    refVec::SVector,
    curr_kmer_freq::Vector{Int},
    consensus_refseq::Seq,
    k::Int64 = 6,
    windowsize::Int64 = 289,
    thr::Union{Int64, Float64} = 33.5,
    buff::Int64 = 50,
    # record_len_thr::Int = 0, # implement in the future: optionally skip record if too short

    mask::UInt64 = unsigned(4095), 
    Nt_bits::DnaBits = NUCLEOTIDE_BITS,
    ScaleFactor::Float64 = 0.166666666666666666666666666666666666666666666667,
    do_align::Bool = true,
    score_model::AffineGapScoreModel{Int64} = AffineGapScoreModel(EDNAFULL, gap_open=-69, gap_extend=-1),
    do_return_dists::Bool = false, 
    do_return_align::Bool = false,
    dist_vec = Float64[], result_align_vec = [],
    resultVec::Vector{FASTA.Record} = FASTA.Record[])

    #get_hit_loci::Bool = false,
    #hit_loci_vec = Int[]
    #genome_pos::Int = 0)

    # get sequence
    seq::Seq = getSeq(record)
    sequence_length::Int = FASTX.FASTA.seqsize(record) # check if using 0.25*sizeof(seq) is better

    # handle length: unfortunately sometimes may skip over sequences in edge cases but should be fine most of the time
    if sequence_length < windowsize; return end 

    #initial operations for the first window  
    fill!(curr_kmer_freq, 0)
    kmer_count!(str = view(seq, 1:windowsize), k = k,
        bins = curr_kmer_freq, mask = mask, Nt_bits = Nt_bits)

    kmerDist::Float64 = ScaleFactor * 0.5 *
        Distances.sqeuclidean(refVec, curr_kmer_freq)

    left_kmer = unsigned(0); for nt in view(seq, 1:k-1)
        left_kmer = (left_kmer << 2) | Nt_bits[nt]
    end

    right_kmer = unsigned(0); for nt in view(seq, windowsize-k+2:windowsize) 
        right_kmer = (right_kmer << 2) | Nt_bits[nt]
    end

    CMI, stop, currminim, i, goal_ind = 2, true, kmerDist, 0, 0
    
    for nt in view(seq, k:sequence_length-windowsize+1); i += 1
        
        @inbounds left_kmer::UInt = ((left_kmer << 2) & mask) | Nt_bits[nt]
        left_ind = left_kmer + 1

        @inbounds right_kmer::UInt = ((right_kmer << 2) & mask) | Nt_bits[seq[i+windowsize]] 
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
                CMI = i
                stop = false
            end
        
        # hit processing 
        elseif !stop
            if CMI > goal_ind; goal_ind = CMI + windowsize - 1
                # alignment and matched unitrange
                left_ind, right_ind = max(CMI-buff,1), min(CMI+windowsize-1+buff,sequence_length)

                if do_align 
                    aligned_obj = pairalign(SemiGlobalAlignment(),
                        view(consensus_refseq, 1:windowsize),
                        view(seq,left_ind:right_ind),score_model)

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
                        " | MatchPos = "*
                        #" | GenomePos = "*
                        " | Len = "*string(last(seq_UnitRange)-first(seq_UnitRange)),
                    view(seq, seq_UnitRange)
                )

                #if get_hit_loci; push!(hit_loci_vec, first(seq_UnitRange)+genome_pos) end 
                push!(resultVec, rec)
                currminim = kmerDist
                stop = true
            end
        end
    end
    #genome_pos += sequence_length
end

function threaded_KmerGMA!(;
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
    score_model::AffineGapScoreModel{Int64} = AffineGapScoreModel(EDNAFULL, gap_open=-69, gap_extend=-1),
    do_return_dists::Bool = false, 
    do_return_align::Bool = false,
    dist_vec = Float64[], result_align_vec = [], 
    resultVec::Vector{FASTA.Record} = FASTA.Record[])
    # get_hit_loci::Bool = false,
    # hit_loci_vec = Int[],
    # genome_pos::Int = 0, 

    # convert refVec to static vector for speed improvement
    refVec = SVector{2 << ((2*k)-1)}(refVec) # takes ab 6 ms
    curr_kmer_freq = zeros(Int, 2 << ((2*k)-1)) # seems like MVector is surprisingly slower

    open(FASTX.FASTA.Reader, genome_path) do reader 
        for record in reader
            record_KmerGMA!( # so the clean approach is to have an input class, but they are slow AF 
                record = record, refVec = refVec, curr_kmer_freq = curr_kmer_freq,
                consensus_refseq = consensus_refseq, # does this create a copy? need to test if theres any improvements if i used subseq
                k = k, windowsize = windowsize, thr = thr, buff = buff,
                mask = mask, Nt_bits = Nt_bits, ScaleFactor = ScaleFactor,
                do_align = do_align, score_model = score_model,
                do_return_dists = do_return_dists, do_return_align = do_return_align,
                dist_vec = dist_vec,
                result_align_vec = result_align_vec, # hit_loci_vec = hit_loci_vec, get_hit_loci = get_hit_loci, genome_pos = genome_pos)
                resultVec = resultVec 
            )
        end
    end
end

export threaded_KmerGMA!

# should use profiler for optimization
# got rid of genome_pos for multithreaded ver bc its a bit complicated and likely slower. the fasta_id_to_cum_len should do the trick.
"""