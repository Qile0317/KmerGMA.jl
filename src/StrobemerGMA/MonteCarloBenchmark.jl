# this is just for the distance not estimation - unfinished, and not generalized for all strobemers
function test_strobe_2_mer_dists(Rkv::String, num_seeds::Int64 = 42, stepsize::Float64 = 0.0125;
    s::Int = 2, w_min::Int = 3, w_max::Int = 5, q::Int = 5)
    reader = open(FASTX.FASTA.Reader, Rkv)
    RKV = strobe_2_mer_count(getSeq(first(reader)), s, w_min, w_max, q)
    close(reader)

    scale_fac = 1/(2*getK(RKV))

    all_dists = Vector{Float64}[]
    for seed in 1:num_seeds
        Random.seed!(seed)
        dists = Float64[]
        for mutation_rate in 0:stepsize:1
            new_vgene = mutate_seq(vgene, mutation_rate)
            push!(dists, scale_fac*Distances.sqeuclidean(
                RKV, strobe_2_mer_count(new_vgene, s, w_min, w_max, q)))
            
        end
        push!(all_dists, dists)
    end
    return all_dists
end

mutation_plot(test_strobe_2_mer_dists("test/Alp_V_ref.fasta", 200; s = 2, w_min = 4, w_max = 10), alpha = 0.05, color = "orange")