using FASTX, Distances, Random
import Base

"""
Function that takes the path of a fasta file / vector of fasta records
and kmer length `k`. Returns a tuple of length 2, where the first value is the reference kmer frequency distribution,
and the second value is the average sequence length.
"""
function gen_ref_ws_cons(reference_seqs, k::Int; get_maxlen = false, Nt_bits::DnaBits = NUCLEOTIDE_BITS)

    is_str = false
    if typeof(reference_seqs) == String
        is_str = true
        reader = open(FASTA.Reader, reference_seqs)
    elseif typeof(reference_seqs) == Vector{FASTX.FASTA.Record}
        reader = reference_seqs # Ref
    else
        stop("invalid input type")
    end

    len, answer, mask = 0, zeros(Float64, 4^k), unsigned(4^k - 1)
    cumulative_nts, maxlen = 0, 0
    curr_profile = Profile(1)

    for record in reader
        len += 1
        curr_len = FASTA.seqsize(record)
        cumulative_nts += curr_len
        if get_maxlen; maxlen = max(maxlen, curr_len) end

        seq = getSeq(record)
        kmer_count!(str=seq,k=k,bins=answer,mask=mask,Nt_bits=Nt_bits)

        lengthen!(curr_profile, curr_len)
        add_consensus!(curr_profile, seq)
    end

    if is_str; close(reader) end
    len = 1/len

    if get_maxlen 
        return answer.*len, Int(round(cumulative_nts*len)), consensus_seq(curr_profile), maxlen
    end
    return answer.*len, Int(round(cumulative_nts*len)), consensus_seq(curr_profile)
end

export gen_ref_ws_cons

# version for RV clustering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# helper function to determine if num is in which range
function get_cluster_index(inp, cutoffs::Vector)
    answer = 1
    for num in cutoffs
        if inp <= num; return answer end
        answer += 1
    end
    return answer
end

export get_cluster_index

"""
    cluster_ref_API(
        reference_seqs::String,
        k::Int;
        maxlen::Int = 0,
        cutoffs::Vector = [7,12,20,25],
        get_dists::Bool = false,
        average_KFV::Vector{Float64} = Float64[],
        Nt_bits::DnaBits = NUCLEOTIDE_BITS)

Function that takes the path of a fasta file / vector of fasta records
and kmer length `k`. Clusters references based on cutoffs, gets the consensus sequences,
makes the new KFVs for each cluster, gets the new windowsizes, 
"""
function cluster_ref_API(reference_seqs::String, k::Int;
    maxlen::Int = 0,
    cutoffs::Vector = [7,12,20,25],
    get_dists::Bool = false,
    average_KFV::Vector{Float64} = Float64[],
    include_avg::Bool = true, # include the overall
    Nt_bits::DnaBits = NUCLEOTIDE_BITS)

    if (maxlen == 0) || (average_KFV == Float64[])
        average_KFV, average_len, average_cons, maxlen = gen_ref_ws_cons(reference_seqs, k; get_maxlen = true, Nt_bits = Nt_bits)
    end

    num_cutoffs = length(cutoffs) + 1

    invalid_vec = [false for _ in 1:num_cutoffs] # false means it IS valid
    lens = [0 for _ in 1:num_cutoffs]

    KFVs, windowsizes = [zeros(4^k) for _ in 1:num_cutoffs], zeros(Int, num_cutoffs)
    consensus_dat_vec::Vector{Any} = Any[Profile(maxlen) for _ in 1:num_cutoffs]
    if get_dists; dists = [] end

    # iterate through records again to get the distance of the reference to the avgKFV
    open(FASTA.Reader, reference_seqs) do reader
        for record in reader
            seq = getSeq(record)
            curr_kmer_dist = kmer_dist(seq, average_KFV, k)
            cluster_ind = get_cluster_index(curr_kmer_dist, cutoffs)
            if get_dists; push!(dists, curr_kmer_dist) end
            
            add_consensus!(consensus_dat_vec[cluster_ind], seq) # could be done during kmer_dist but ehh...
            windowsizes[cluster_ind] += length(seq)
            lens[cluster_ind] += 1
            kmer_count!(str = seq, k = k, bins = KFVs[cluster_ind], mask = unsigned((4^k)-1), Nt_bits = Nt_bits)
        end
    end

    consensus_vec = [dna"" for _ in 1:num_cutoffs]
    for i in 1:num_cutoffs
        if lens[i] != 0
            KFVs[i] ./= lens[i]
            windowsizes[i] = Int(round(windowsizes[i]/lens[i]))
            consensus_vec[i] = consensus_seq(consensus_dat_vec[i])[1:windowsizes[i]]
        else
            invalid_vec[i] = true
        end
    end

    if include_avg 
        push!(invalid_vec, false)
        push!(KFVs, average_KFV)
        push!(windowsizes, average_len)
        push!(consensus_vec, average_cons)
    end

    if get_dists
        return KFVs, windowsizes, consensus_vec, invalid_vec, dists 
    end
    return KFVs, windowsizes, consensus_vec, invalid_vec
end

export cluster_ref_API

function eliminate_null_params()