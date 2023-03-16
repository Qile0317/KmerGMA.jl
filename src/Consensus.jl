using BioSequences
import Base

const Num_dict = Dict{Int, DNA}(1 => DNA_A, 2 => DNA_C, 3 => DNA_G, 4 => DNA_T)

mutable struct Profile
    vecs::Vector{Vector{Int64}}
    len::Int64

    Profile(len::Int64) = new([zeros(Int64, len) for _ in 1:4], len)
    Base.getindex(p::Profile, nt::DNA) = p.vecs[NUCLEOTIDE_BITS[nt]+1]
end

export Profile

function add_consensus!(vec::Profile, seq)
    for i in eachindex(seq) 
        vec[seq[i]][i] += 1 
    end
end

export add_consensus!

function lengthen!(p::Profile, new_len::Int64)
    if new_len > p.len
        for i in 1:4
            for _ in 1:(new_len - p.len)
                push!(p.vecs[i],0)
            end
        end
        p.len = new_len
    end
end

export lengthen!

function consensus_seq(consensus::Profile, )
    curr_max, curr_cons_max = dna"A"^(consensus.len), copy(consensus.vecs[1])
    for i in 2:4
        for j in 1:consensus.len
            if consensus.vecs[i][j] > curr_cons_max[j]
                curr_cons_max[j] = consensus.vecs[i][j]
                curr_max[j] = Num_dict[i]
            end
        end
    end
    return curr_max
end

export consensus_seq