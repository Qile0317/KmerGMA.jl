using FASTX

export findGenes
export findGenes_cluster_mode
export write_results

# helper function for warnings
function warn_helper(k::Int, do_return_dists::Bool)
    if k < 5; @warn "Such a low k value of $k likely won't yield the most accurate results" end
    if do_return_dists; @warn "Setting do_return_dists to true may be very memory intensive" end 
end

"""
    KmerGMA.findGenes(;
        genome_path::String,
        ref_path::String,
        k::Int = 6,
        KmerDistThr::Union{Int64, Float64} = 0,
        buffer::Int64 = 50,
        do_align::Bool = true,
        gap_open_score::Int = -69,
        gap_extend_score::Int = -1,
        do_return_dists::Bool = false,
        do_return_hit_loci::Bool = false,
        do_return_align::Bool = false
        verbose::Bool = true)

The main API to conduct homology searching in a genome, using a kmer-based sequence similarity metric, against a a reference sequence set. For example, the set of all germline V genes of a mammal.
Returns the approximate matches as FASTA record vector WITHIN A VECTOR of length 1 in the default configuration of the parameters. The descriptions of the record contain information about the match. 
The format of the description appears as so:

   "Identifier | dist = a | MatchPos = b:c | GenomePos = e | Len = f"

Where `Identifier` is the contig ID of the genome where the current hit was found. `dist` is the kmer similarity, and `MatchPos` is the unitrange for the match in the contig.
`GenomePos` is the cumulative nucleotides that have been iterated over until the current contig, and `Len` is simply the length of the hit length.

...
# Arguments
- `genome_path`: a string that should be the path of the genome to conduct homology searching on.
- `ref_path`: the file location of a fasta file containing the reference sequence(s). The references should be very similar in length.

# Optional Arguments
- `k::Int64 = 6`: the kmer length to use for approximate matching. Generally it should probably remain between 5 - 10 and probably does not have a profound impact on how the matches are found. Check out the pre-print for more information
- `KmerDistThr::Int64 = 0`: the Kmer Distance distance threshold for sequence matches to the query reference sequence set. Out of the context of the algorithm, lower values mean matches have to be more similar to the references. If left as 0, it is automatically computed. Once again the pre-print has information on this argument.
- `buffer::Int64 = 50`: the amount of nucleotides left and right of a matched sequence that should be added to the returned fasta sequences, as KmerGMA is a heuristic. If `do_align` is set to `true` the buffer will be included for a more accurate alignment.
- `verbose::Bool = true` Indicates whether to show info in the REPL about the the progress of the processing
- `do_align`: Whether to align the hits+buffer region to the consensus sequence of the references. Highly recommended to keep as `true`
- `gap_open_score::Int = -69`: The `gap_open` scoring paramater for the `BioAlignments` semi global pairwise aligner function. The number should ideally be kept low to heavily penalize gap extensions for most use-cases.
- `gap_extend_score::Int = -1`: The `gap_extend` scoring parameter for the `BioAlignments` semi global pairwise aligner function. Can be kept to the default `-1` as long as the `gap_open` score is low.
- `do_return_dists`: boolean to indicate whether the kmer distances along every window along the genome should be returned in a vector. (intensive memory consumption when genomes are large)
- `do_return_hit_loci`: if `true`, will return an additional vector of the position within the genomic sequences of each hit, corresponding to the index in the hit vector.
- `do_return_align`: if `true`, will return an additional vector of alignment object of each hit to the consensus reference sequence.
- `KmerDist_threshold_buffer::Real = 8.0`: a value to determine approximately how much kmer distance a hit should be lower than any random non-hit sequence. Keep in mind that kmerdistance approximates edit distance for mutations better than indels.
...

The last three arguments would add term to the output. The output vector would incorporate the respective vectors in the same order of priority if any of the parameters are true.

Note: Playing with the `KmerDistThr` argument could return more or less matches every time. 
"""
function findGenes(; genome_path::String, ref_path::String,
    k::Int = 6, KmerDistThr::Union{Int64, Float64} = 0, buffer::Int64 = 50,
    do_align::Bool = true, gap_open_score::Int = -69, gap_extend_score::Int = -1,
    do_return_dists::Bool = false, do_return_hit_loci::Bool = false, do_return_align::Bool = false,
    verbose::Bool = true, KmerDist_threshold_buffer::Real = 8.0)

    if verbose; @info "pre-processing references and parameters..." end 
    warn_helper(k, do_return_dists)

    RV, windowsize, consensus_refseq = gen_ref_ws_cons(ref_path, k)
    if k >= windowsize; error("the average reference sequence length $windowsize exceeds/is equal to the chosen kmer length $k. please reduce k. ") end

    estimated_optimal_KmerDistThr = estimate_optimal_threshold(RV, windowsize, buffer = KmerDist_threshold_buffer)
    if KmerDistThr == 0
        KmerDistThr = estimated_optimal_KmerDistThr
    elseif KmerDistThr < estimated_optimal_KmerDistThr
        @warn "The kmer distance threshold $KmerDistThr for k = $k is likely too high, and can result in many false positives"
    end
    
    hit_vector = FASTX.FASTA.Record[]
    dist_vec, hit_loci_vec, alignment_vec = Float64[], Int[], []  
    cumulative_length_in_genome = 0

    if verbose; @info "initializing iteration..." end 
    ac_gma_testing!(
        genome_path = genome_path, refVec = RV, consensus_refseq = consensus_refseq,
        k = k, windowsize = windowsize, thr = KmerDistThr, buff = buffer,
        mask = unsigned(4^k - 1), ScaleFactor = 1/k, 
        do_align = do_align, gap_open_score = gap_open_score, gap_extend_score = gap_extend_score,
        do_return_dists = do_return_dists,
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

"""
    KmerGMA.findGenes_cluster_mode(;
        genome_path::String,
        ref_path::String,
        cluster_cutoffs = [7,12,20,25],
        k::Int = 6, 
        KmerDistThrs = Float64[0], 
        buff::Int64 = 50,
        do_align::Bool = true,
        gap_open_score::Int = -200,
        gap_extend_score::Int = -1,
        do_return_align::Bool = false,
        do_return_dists::Bool = false,
        do_return_hit_loci::Bool = false,
        verbose::Bool = true,
        kmerDist_threshold_buffer::Real = 7
        )

A slower (`O(mn)`) alternative to `KmerGMA.findGenes` to conduct homology searching in a genome, against a a reference sequence set, using a kmer-based sequence similarity metric.
For example, the set of all germline V genes of a mammal.

Returns the approximate matches as FASTA record vector WITHIN A VECTOR of length 1 in the default configuration of the parameters. The descriptions of the record contain information about the match. 
The format of the description appears as so:

    "Identifier | dist = a | KFV = b | MatchPos = c:d | GenomePos = e | Len = f"

Where `Identifier` is the contig ID of the genome where the current hit was found. `dist` is the kmer similarity,
`KFV` is the reference kmer frequency vector that a hit was matched against (more info below), `MatchPos` is the unitrange for the match in the contig.
`GenomePos` is the cumulative nucleotides that have been iterated over until the current contig, and `Len` is simply the length of the hit sequence

...
# Arguments
- `genome_path`: a string that should be the path of the genome to conduct homology searching on.
- `ref_path`: the file location of a fasta file containing the reference sequence(s). The references should be very similar in length.

# Optional Arguments
- `cluster_cutoffs = [7,12,20,25]`: Important determination of the cutoff points to construct subsets of similar reference sequences. The cutoff points refer to the kmer distance against the reference. Usually, the more cutoffs may be better but the default cutoffs have been determined to be relatively robust.
- `k::Int64 = 6`: the kmer length to use for approximate matching. Generally it should probably remain between 5 - 10 and probably does not have a profound impact on how the matches are found. Check out the pre-print for more information
- `KmerDistThrs = [0]`: the Kmer Distance distance thresholds for each sequence matche to the query reference sequence set. Out of the context of the algorithm, lower values mean matches have to be more similar to the references. If left as 0, it is automatically computed. Once again a pre-print has information on this argument.
- `buffer::Int64 = 100`: the amount of nucleotides left and right of a matched sequence that should be added to the returned fasta sequences, as KmerGMA is a heuristic. If `do_align` is set to `true` the buffer will be included for a more accurate alignment.
- `do_align`: Whether to align the hits+bufer region to the consensus sequence of the references. Highly recommended to keep as `true`
- `gap_open_score::Int = -200`: The `gap_open` scoring paramater for the `BioAlignments` semi global pairwise aligner function. The number should ideally be kept low to heavily penalize gap extensions for most use-cases.
- `gap_extend_score::Int = -1`: The `gap_extend` scoring parameter for the `BioAlignments` semi global pairwise aligner function. Can be kept to the default `-1` as long as the `gap_open` score is low.
- `do_return_align`: if `true`, will return an additional vector of alignment object of each hit to the consensus reference sequence.
- `do_return_dists`: boolean to indicate whether the kmer distances along every window along the genome should be returned in a vector. (intensive memory consumption when genomes are large)
- `do_return_hit_loci`: if `true`, will return an additional vector of the position within the genomic sequences of each hit, corresponding to the index in the hit vector.
- `verbose::Bool = true` Indicates whether to show info in the REPL about the the progress of the processing
- `kmerDist_threshold_buffer::Real = 7`: a value to determine approximately how much kmer distance a hit should be lower than any random non-hit sequence. Keep in mind that kmerdistance approximates edit distance for mutations better than indels. Additionally, the false positive rate would likely increase as the value goes too low. Unfortunately this value is currently only optimized for `k = 6`. something around `20` seem to work well for `k = 5`. 
...

The last three arguments would add terms to the output. When `verbose` is true the exact outputs are even stated. 
The output vector would incorporate the respective vectors in the same order of priority if any of the parameters are true.

There are some more details to be added in future releases.
"""
function findGenes_cluster_mode(; genome_path::String, ref_path::String,
    cluster_cutoffs = [7,12,20,25], # cluster cutoff visualization function should be made
    k::Int = 6, KmerDistThrs::Kfv = Float64[0.0], buffer::Int64 = 100,
    
    do_align::Bool = true, gap_open_score::Int = -200, gap_extend_score::Int = -1,

    do_return_dists::Bool = false, do_return_hit_loci::Bool = false, do_return_align::Bool = false,
    verbose::Bool = true,
    kmerDist_threshold_buffer::Real = 7) # (It is dependent on the average variance of of kmer distance so a function can do that)

    if verbose; @info "pre-processing references and parameters..." end 
    warn_helper(k, do_return_dists)

    RVs, windowsizes, consensus_refseqs, invalids = cluster_ref_API(ref_path, k; cutoffs = cluster_cutoffs)
    RVs, windowsizes, consensus_refseqs = eliminate_null_params(RVs, windowsizes, consensus_refseqs, invalids); invalids = nothing

    if k >= minimum(windowsizes); error("some/all of the average reference sequence lengths exceeds/is equal to the chosen kmer length $k. please reduce k. ") end

    # distance thresholds
    estimated_optimal_KmerDistThrs = estimate_optimal_threshold(RVs, windowsizes; buffer = kmerDist_threshold_buffer)
    if first(KmerDistThrs) == 0 # user indicates they want it to be automated
        KmerDistThrs = estimated_optimal_KmerDistThrs
    else
        num_warn_str, ind_warn_str = "", ""
        for (i, num) in enumerate(KmerDistThrs)
            if num > estimated_optimal_KmerDistThrs[i]
                ind_warn_str *= "$i, "
                num_warn_str *= "$num, "
            end
        end
        if ind_warn_str != ""
            @warn "The kmer distance thresholds $KmerDistThrs"*" at index/indicies $ind_warn_str"[1:end-2]*" for k = $k is potentially too high, and may result in more false positives."
        end
    end; estimated_optimal_KmerDistThrs = nothing # hopefully saves teeny bit of memory
    
    hit_vector = FASTX.FASTA.Record[]
    hit_loci_vec, alignment_vec = Int[], SubAlignResult[]  
    dist_vec_vec = [Float64[] for _ in 1:length(windowsizes)]

    if verbose; @info "initializing iteration..." end 
    Omn_KmerGMA!(genome_path = genome_path, refVecs = RVs,
        windowsizes = windowsizes, consensus_seqs = consensus_refseqs,
        resultVec = hit_vector, k = k,
        ScaleFactor = 1/k, mask = unsigned(4^k -1),
        thr_vec = KmerDistThrs,
        buff = buffer,

        align_hits = do_align,
        gap_open_score = gap_open_score,
        gap_extend_score = gap_extend_score,
        get_aligns = do_return_align,
        
        get_hit_loci = do_return_hit_loci, 
        hit_loci_vec = hit_loci_vec,
        align_vec = alignment_vec,
        do_return_dists = do_return_dists, dist_vec_vec = dist_vec_vec)

    info_str = "genome mining completed successfully, returning vector of: vector of hits"
    output_vector = Any[hit_vector]
    if do_return_hit_loci; push!(output_vector, hit_loci_vec); info_str *= ", vector of hit locations" end
    if do_return_align; push!(output_vector, alignment_vec); info_str *= ", vector of alignments" end 
    if do_return_dists; push!(output_vector, dist_vec_vec) ; info_str *= ", vector of vectors of kmer distances along the genome"end
    
    if verbose; @info info_str; @info "To write the results, use `KmerGMA.write_results`" end
    return output_vector
end

"""
    write_results(KmerGMA_result_vec::Vector{FASTX.FASTA.Record}, file_path::String, width::Int64 = 95)

Writes all FASTA Records (from the FASTX package) from a vector into a fasta file location `file_path`.
`width` is an optional argument that indicates the maximum width of sequences written to the file per line.
"""
function write_results(KmerGMA_result_vec::Vector{FASTX.FASTA.Record}, file_path::String, width::Int64 = 95)
    FASTA.Writer(open(file_path, "a"), width = width) do writer
        for hit in KmerGMA_result_vec
            write(writer, hit)
        end
    end
    @info "writing complete"
end