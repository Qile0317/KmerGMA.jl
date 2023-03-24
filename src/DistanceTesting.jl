using BioSequences, Random

export estimate_optimal_threshold

# looks at average distance of random sequence to the RV and returns value just below
function estimate_optimal_threshold(RV::Kfv, average_length::Int64;
    seed::Int64 = 42, num_trials::Int64 = 100, buffer::Real = 8)

    Random.seed!(seed); k = Int(log(4, length(RV)))
    cumulative_dist_to_RV = 0
    for trial in 1:num_trials
        cumulative_dist_to_RV += kmer_dist(randdnaseq(average_length), RV, k)
    end
    return (cumulative_dist_to_RV/num_trials) - buffer
end

function estimate_optimal_threshold(RV_vec::Vector{Vector{Float64}}, average_length_vec::Vector{Int64};
    seed::Int64 = 42, num_trials::Int64 = 100, buffer::Real = 8)

    Random.seed!(seed); k = Int(log(4, length(first(RV_vec))))
    thresholds = Float64[]
    for i in eachindex(RV_vec)
        RV, cumulative_dist_to_RV = RV_vec[i], 0
        for trial in 1:num_trials
            cumulative_dist_to_RV += kmer_dist(randdnaseq(average_length_vec[i]), RV, k)
        end
        push!(thresholds, (cumulative_dist_to_RV/num_trials) - buffer)
    end
    return thresholds
end

# to have a better buffer, the variance should be determined and take the lower bound - some additional val.
# more extensive kmer distance testing needs to be done
# additionally the RV can be mutated

const mutation_dict = Dict{DNA, Vector{DNA}}(
    DNA_A => DNA[DNA_C, DNA_G, DNA_T],
    DNA_C => DNA[DNA_A, DNA_G, DNA_T],
    DNA_G => DNA[DNA_C, DNA_A, DNA_T],
    DNA_T => DNA[DNA_C, DNA_G, DNA_A]
)

function mutate_seq(seq::Seq, mut_rate::Real)
    newseq = copy(seq)
    for i in 1:length(seq)
        if rand(1)[1] <= mut_rate
            newseq[i] = rand(mutation_dict[seq[i]])
        end
    end
    return newseq
end

function mutate_seq!(seq::Seq, mut_rate::Real)
    for i in 1:length(seq)
        if rand(1)[1] <= mut_rate
            seq[i] = rand(mutation_dict[seq[i]])
        end
    end
end

# plot substitution testing results
function mutation_plot(test_res::Vector{Vector{Float64}}; stepsize::Float64 = 0.0125, alpha::Float64 = 0.05)
    x_axis_vec = [x for x in 0:0.0125:1]
    return scatter([x_axis_vec for _ in 1:length(test_res)], test_res,
        alpha = alpha, color =:black, label = nothing,
        xlabel = "mutation rate", ylabel = "distance")
end