using BioSequences, FASTX, Distances, BioAlignments

const InputConsts = NamedTuple{(:genome_path, :refVec, :consensus_refseq, :k, :windowsize, :thr, :buff, :mask, :Nt_bits, :ScaleFactor, :do_align, :score_model, :do_return_dists, :do_return_align, :get_hit_loci),
    Tuple{String, Vector{Float64}, LongSequence{DNAAlphabet{4}}, Int64, Int64, Float64, Int64, UInt64, Dict{DNA, UInt64}, Float64, Bool, AffineGapScoreModel{Int64}, Bool, Bool, Bool}}

function init_InputConsts(;
    genome_path::String,
    refVec::Vector{Float64},
    consensus_refseq::Seq,
    k::Int64 = 6,
    windowsize::Int64 = 289,
    thr::Union{Int64, Float64} = 33.5,
    buff::Int64 = 50,
    mask::UInt64 = unsigned(4095), 
    Nt_bits::DnaBits = NUCLEOTIDE_BITS,
    ScaleFactor::Float64 = 0.0833333333333333333333,
    do_align::Bool = true,
    score_model::AffineGapScoreModel{Int64} = AffineGapScoreModel(EDNAFULL, gap_open=-5, gap_extend=-1),
    do_return_dists::Bool = false, # dangerous for memory if true
    do_return_align::Bool = false,
    get_hit_loci::Bool = false
    )
    return InputConsts((genome_path,refVec,consensus_refseq,
        k,windowsize,thr,buff,mask,Nt_bits,ScaleFactor,do_align,
        score_model,do_return_dists,do_return_align,get_hit_loci))
end

export init_InputConsts

# sooo... from benchmarktools it seems using the namedTuple actually worsened the performance...
# additionally, lots of memory is consumed because the records become stored in memory

function ac_gma_testing!(; inp::InputConsts,
    curr_kmer_freq::Vector{Float64} = zeros(4096),
    dist_vec = Float64[], result_align_vec = [], hit_loci_vec = Int[],
    genome_pos::Int = 0, resultVec::Vector{FASTA.Record} = FASTA.Record[]
    )
    open(FASTX.FASTA.Reader, inp.genome_path) do reader 
        for record in reader 
            k = inp.k
            seq::Seq = getSeq(record)
            proceed::Bool, goal_ind::Int64= true, 0

            #edge case
            sequence_length = length(seq)
            if sequence_length < inp.windowsize; return end

            #initial operations for the first window  
            fill!(curr_kmer_freq, 0)
            kmer_count!(str = view(seq, 1:inp.windowsize), k = k,
                bins = curr_kmer_freq, mask = inp.mask, Nt_bits = inp.Nt_bits)
            currSqrEuc = Distances.sqeuclidean(inp.refVec, curr_kmer_freq)

            #initializing variables
            CMI, stop, currminim = 2, true, currSqrEuc

            left_kmer = unsigned(0)
            for c in seq[1:k-1]
                left_kmer = (left_kmer << 2) + inp.Nt_bits[c]
            end

            right_kmer = unsigned(0)
            for c in seq[inp.windowsize-k+1:inp.windowsize-1]
                right_kmer = (right_kmer << 2) + inp.Nt_bits[c]
            end

            for i in 1:(sequence_length-inp.windowsize)
                # first kmer
                left_kmer = ((left_kmer << 2) & inp.mask) + inp.Nt_bits[seq[i+k-1]]
                left_ind = -~left_kmer
                @views currSqrEuc -= (inp.refVec[left_ind]-curr_kmer_freq[left_ind])^2 
                curr_kmer_freq[left_ind] -= 1
                @views currSqrEuc += (inp.refVec[left_ind]-curr_kmer_freq[left_ind])^2

                # last kmer + 1bp
                right_kmer = ((right_kmer << 2) & inp.mask) + inp.Nt_bits[seq[i+inp.windowsize-1]]
                right_ind = -~right_kmer
                @views currSqrEuc -= (inp.refVec[right_ind]-curr_kmer_freq[right_ind])^2
                curr_kmer_freq[right_ind] += 1
                @views currSqrEuc += (inp.refVec[right_ind]-curr_kmer_freq[right_ind])^2

                # convert to kmer Distance 
                kmerDist = currSqrEuc * inp.ScaleFactor 
                if inp.do_return_dists; push!(dist_vec, kmerDist) end 

                # minima finder
                if kmerDist < inp.thr
                    if kmerDist < currminim
                        currminim = kmerDist
                        CMI = i
                        stop = false
                    end
                
                # hit processing 
                elseif !stop
                    if !proceed 
                        if CMI > goal_ind
                            proceed = true
                        end
                    else # operations here aren't optimized but for actual long sequences the time is negligible
                        goal_ind = CMI + inp.windowsize - 1
                        proceed = false 
                        
                        # alignment and matched unitrange
                        left_ind, right_ind = max(CMI-inp.buff,1), min(CMI+inp.windowsize-1+inp.buff,sequence_length)

                        if inp.do_align 
                            aligned_obj = pairalign(SemiGlobalAlignment(),inp.consensus_refseq,view(seq,left_ind:right_ind),inp.score_model)
                            if inp.do_return_align; push!(result_align_vec, aligned_obj) end
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
                                " | GenomePos = $genome_pos",
                            view(seq, seq_UnitRange)
                        )

                        if inp.get_hit_loci; push!(hit_loci_vec, first(seq_UnitRange)+genome_pos) end 
                        push!(resultVec, rec)
                        currminim = currSqrEuc
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

"""
a = (genome_path = "hi", refVec = [0.0,0.1], consensus_refseq = dna"ATG",k= 6,windowsize = 289,
           thr = 33.5,
           buff = 50,
           curr_kmer_freq = zeros(4096),
           mask = unsigned(4095),
           Nt_bits = NUCLEOTIDE_BITS,
           ScaleFactor = 0.0833333333333333333333,do_align = true,vscore_model = AffineGapScoreModel(EDNAFULL, gap_open=-5, gap_extend=-1),

           do_return_dists = false, dist_vec = Float64[],do_return_align = false,
           result_align_vec = [],get_hit_loci = false,hit_loci_vec = Int[],genome_pos = 0, resultVec = FASTA.Record[])
"""

# alternative approach for multiple RVs that account for duplicates. Bugged at the moment.

function Omn_gma_testing(;
    genome_path::String,
    genome_pos::Int = 0,
    k::Int64 = 6,
    refVecs::Vector{Vector{Float64}} = [reference_KFV, ref_cl2_KFV,single_ref_KFV],
    windowsize::Int64 = 298, #could be a vec
    thr::Union{Int64, Float64} = 30, # could be a vec
    buff::Int64 = 50, # could be a vec
    curr_kmer_freq::Vector{Float64} = zeros(4096), # no bin for N!!!
    mask::UInt64 = unsigned(4095), 
    Nt_bits::Dict{DNA, UInt64} = NUCLEOTIDE_BITS,
    ScaleFactor::Float64 = 0.0833333333333333333333,
    score_model::AffineGapScoreModel{Int64} = AffineGapScoreModel(EDNAFULL, gap_open=-69, gap_extend=-1),
    cons_gene::LongSequence{DNAAlphabet{4}} = vgene,
    do_align::Bool = true,
    do_return_dists::Bool = false,
    do_return::Bool = true,
    resultVec = FASTA.Record[],
    result_align_vec = [],
)
    
    len_KFVs = length(refVecs)
    currSqrEuc_vec = Float64[0.0 for _ in 1:len_KFVs]
    if do_align; result_align_vec = [] end
    if do_return_dists; dist_vec_vec = [Float64[] for _ in 1:len_KFVs] end 

    CMIs = [1 for _ in 1:len_KFVs]
    stops = [true for _ in 1:len_KFVs]

    # begin iteration (could be made into a function)
    open(FASTX.FASTA.Reader, genome_path) do reader 
        for record in reader 
            
            # Operation on first window 
            seq::LongSequence{DNAAlphabet{4}} = KmerGMA.getSeq(record)

            sequence_length::Int64 = length(seq)
            if sequence_length < windowsize; continue end
            
            goal_ind_range, goal_ind = 0:0, 0

            #initial operations for the first window  
            fill!(curr_kmer_freq, 0)

            KmerGMA.kmer_count!(str = view(seq, 1:windowsize), k = k,
            bins = curr_kmer_freq, mask = mask, Nt_bits = Nt_bits)

            # for each KFV
            for KFV_ind in 1:len_KFVs
                currSqrEuc_vec[KFV_ind] = Distances.sqeuclidean(refVecs[KFV_ind], curr_kmer_freq)
            end

            #initializing variables for distance
            fill!(CMIs, 1); fill!(stops, true)
            curr_mins = [num for num in currSqrEuc_vec]

            left_kmer = unsigned(0)
            for c in seq[1:k-1]
                left_kmer = (left_kmer << 2) + Nt_bits[c]
            end

            right_kmer = unsigned(0)
            for c in seq[windowsize-k+1:windowsize-1]
                right_kmer = (right_kmer << 2) + Nt_bits[c]
            end

            # iteration 
            for i in 1:(sequence_length-windowsize) 

                left_kmer = ((left_kmer << 2) & mask) + Nt_bits[seq[i+k-1]]
                left_ind = -~left_kmer

                right_kmer = ((right_kmer << 2) & mask) + Nt_bits[seq[i+windowsize-1]]
                right_ind = -~right_kmer

                for ind in 1:len_KFVs  
                    currSqrEuc_vec[ind] -= (refVecs[ind][left_ind] - curr_kmer_freq[left_ind])^2 
                    currSqrEuc_vec[ind] += (refVecs[ind][left_ind] - curr_kmer_freq[left_ind] + 1)^2

                    currSqrEuc_vec[ind] -= (refVecs[ind][right_ind] - curr_kmer_freq[right_ind])^2
                    currSqrEuc_vec[ind] += (refVecs[ind][right_ind] - curr_kmer_freq[right_ind] - 1)^2
                end

                curr_kmer_freq[left_ind] -= 1
                curr_kmer_freq[right_ind] += 1

                for ind in 1:len_KFVs
                    kmerDist = currSqrEuc_vec[ind] * ScaleFactor 
                    if do_return_dists; push!(dist_vec_vec[ind], kmerDist) end

                    # minima finder (worse matches may sometimes be used, needs fix)
                    if kmerDist < thr
                        if kmerDist < curr_mins[ind]
                            curr_mins[ind] = kmerDist
                            CMIs[ind] = i
                            stops[ind] = false
                        end
                    elseif !stops[ind] # if kmerDist is above/equal threshold, check if do stop
                        CMI = CMIs[ind]
                        #print(" $CMI")

                        if ((CMI in goal_ind_range) == false)
                            goal_ind_range = CMI-windowsize+1:CMI + windowsize - 1 # is something wrong here?

                            # alignment and matched unitrange
                            left_ind, right_ind = max(CMI-buff,1), min(CMI+windowsize-1+buff,sequence_length)

                            if do_align 
                                aligned_obj = pairalign(SemiGlobalAlignment(), cons_gene, view(seq, left_ind:right_ind), score_model)
                                if do_return; push!(result_align_vec, aligned_obj) end 
                                aligned_UnitRange = cigar_to_UnitRange(aligned_obj)
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
                        end
                        curr_mins[ind] = currSqrEuc_vec[ind]
                        stops[ind] = true 
                    end
                end
            end
            genome_pos += sequence_length
        end
    end
    if do_return 
        if do_align
            if do_return_dists; return resultVec, result_align_vec, dist_vec_vec end
            return resultVec, result_align_vec
        end
        return resultVec
    end
end