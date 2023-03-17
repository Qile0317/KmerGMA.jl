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
    seed::Int64 = 42, num_trials::Int64 = 100, buffer::Real = 10)

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

# unfinished - use to visualize
function mutation_plot(RV::Kfv; seed::Int64 = 42, num_trials::Int64 = 50)
    println("placeholder")
end