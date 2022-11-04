module KmerGMA
    using Plots, BioSequences, FASTX, Distances, ProgressMeter #WebBlast, Statistics, StringDistances, StatsBase, StatsPlots, DataFrames
    include("simpleExplore.jl")
    include("ExactMatch.jl")
    include("RefGen.jl")
    include("GMA.jl")
    include("API.jl")
    include("BLAST.jl")
end
