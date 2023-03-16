using FASTX

# version without clustering

"""
    KmerGMA.findGenes(;
        genome_path::String,
        ref_path::String,
        do_cluster::Bool = false, # to be implemented
        k::Int = 6,
        KmerDistThr::Union{Int64, Float64} = 0,
        buffer::Int64 = 50,
        do_align::Bool = true,
        do_return_dists::Bool = false,
        do_return_hit_loci::Bool = false,
        do_return_align::Bool = false)

The main API to conduct homology searching in a genome, using a kmer-based sequence similarity metric, against a a reference sequence set. For example, the set of all germline V genes of a mammal.
Returns the approximate matches as FASTA record vector WITHIN A VECTOR of length 1 in the default configuration of the parameters. The descriptions of the record contain information about the match. 
The format of the description appears as so:

   "Identifier | dist = a | MatchPos = b:c | GenomePos = e"

Where `Identifier` is the contig ID of the genome where the current hit was found. `dist` is the kmer similarity, and `MatchPos` is the unitrange for the match in the contig.
`GenomePos` is the cumulative nucleotides that have been iterated over until the current contig.

...
# Arguments
- `genome_path`: a string that should be the path of the genome to conduct homology searching on.
- `ref_path`: the file location of a fasta file containing the reference sequence(s). The references should be very similar in length.

# Optional Arguments
- `do_cluster::Bool = false`: work in progress, has no functionality at the moment
- `k::Int64 = 6`: the kmer length to use for approximate matching. Generally it should probably remain between 5 - 10 and probably does not have a profound impact on how the matches are found. Check out the pre-print for more information
- `KmerDistThr::Int64 = 0`: the Kmer Distance distance threshold for sequence matches to the query reference sequence set. Out of the context of the algorithm, lower values mean matches have to be more similar to the references. If left as 0, it is automatically computed. Once again the pre-print has information on this argument.
- `buffer::Int64 = 50`: the amount of nucleotides left and right of a matched sequence that should be added to the returned fasta sequences, as KmerGMA is a heuristic
- `do_align`: Whether to align the hits+bufer region to the consensus sequence of the references. Highly recommended to keep as `true`
- `do_return_dists`: dangerous boolean to indicate whether the kmer distances along every window alon ghte genome should be returned in a vector
- `do_return_hit_loci`: if `true`, will return an additional vector of the position within the genomic sequences of each hit, corresponding to the index in the hit vector.
- `do_return_align`: if `true`, will return an additional vector of alignment object of each hit to the consensus reference sequence.
...

The last three arguments would add term to the output. The output vector would incorporate the respective vectors in the same order of priority if any of the parameters are true.

Note: Playing with the `KmerDistThr` argument could return more or less matches every time. 
"""
function findGenes(; genome_path::String, ref_path::String, do_cluster::Bool = false, # to be implemented
    k::Int = 6, KmerDistThr::Union{Int64, Float64} = 0, buffer::Int64 = 50, do_align::Bool = true,
    do_return_dists::Bool = false, do_return_hit_loci::Bool = false, do_return_align::Bool = false)

    if k < 5; @warn "Such a low k value of $k likely won't yield the most accurate results" end
    if do_return_dists; @warn "Setting do_return_dists to true may be very memory intensive" end 

    RV, windowsize, consensus_refseq = gen_ref_ws_cons(ref_path, k)

    estimated_optimal_KmerDistThr = estimate_optimal_threshold(RV, windowsize)
    if KmerDistThr == 0
        KmerDistThr = estimated_optimal_KmerDistThr
    elseif KmerDistThr < estimated_optimal_KmerDistThr
        @warn "The kmer distance threshold $KmerDistThr for k = $k is likely too high, and can result in many false positives"
    end
    
    KmerGMA_constant_parameters = init_InputConsts(
        genome_path = genome_path, refVec = RV, consensus_refseq = consensus_refseq,
        k = k, windowsize = windowsize, thr = KmerDistThr, buff = buffer,
        mask = unsigned(4^k - 1), ScaleFactor = 1/(2*k), 
        do_align = do_align, do_return_dists = do_return_dists,
        do_return_align = do_return_align, get_hit_loci = do_return_hit_loci);
    
    hit_vector = FASTX.FASTA.Record[]
    dist_vec, hit_loci_vec, alignment_vec = Float64[], Int[], []  
    curr_KFV, cumulative_length_in_genome = zeros(4^k), 0

    ac_gma_testing!(inp = KmerGMA_constant_parameters,
        curr_kmer_freq = curr_KFV, dist_vec = dist_vec,
        result_align_vec = alignment_vec,
        hit_loci_vec = hit_loci_vec, genome_pos = cumulative_length_in_genome,
        resultVec = hit_vector)

    # return results
    output_vector = Any[hit_vector]
    if do_return_hit_loci; push!(output_vector, hit_loci_vec) end
    if do_return_align; push!(output_vector, alignment_vec) end 
    if do_return_dists; push!(output_vector, dist_vec) end
    return output_vector
end

export findGenes

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
    println("complete")
end

export write_results