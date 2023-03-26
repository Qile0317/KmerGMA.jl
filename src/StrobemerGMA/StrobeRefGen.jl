# this function is also in ReferenceGeneration.jl !

# w_min = 4 and w_max = 7 might work a bit better intuitively
function gen_ref_ws_cons(reference_seqs;
    s::Int = 2, w_min::Int = 3, w_max::Int = 5, q::Int = 5,
    get_maxlen = false)

    is_str = false
    if typeof(reference_seqs) == String
        is_str = true
        reader = open(FASTA.Reader, reference_seqs)
    elseif typeof(reference_seqs) == Vector{FASTX.FASTA.Record}
        reader = reference_seqs # Ref
    else
        stop("invalid input type")
    end

    len, answer = 0, zeros(Float64, 2 << ((4*s)-1))
    cumulative_nts, maxlen = 0, 0
    curr_profile = Profile(1)

    for record in reader
        len += 1
        curr_len = FASTA.seqsize(record)
        cumulative_nts += curr_len
        if get_maxlen; maxlen = max(maxlen, curr_len) end

        seq = getSeq(record)
        ungapped_strobe_2_mer_count!(seq, answer,
            s, w_min, w_max, q)

        lengthen!(curr_profile, curr_len)
        add_consensus!(curr_profile, seq)
    end

    if is_str; close(reader) end
    len = 1/len

    if get_maxlen 
        return SVector{2 << ((4*s)-1)}(answer.*len), Int(round(cumulative_nts*len)), consensus_seq(curr_profile), maxlen
    end
    return SVector{2 << ((4*s)-1)}(answer.*len), Int(round(cumulative_nts*len)), consensus_seq(curr_profile)
end