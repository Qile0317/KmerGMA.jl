using BioSequences, Random

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

export estimate_optimal_threshold

# to have a better buffer, the variance should be determined and take the lower bound - some additinoal val.
# will put in the mutation plotting stuff here

function mutation_plot(RV::Kfv; seed::Int64 = 42, num_trials::Int64 = 50)
    println("placeholder")
end