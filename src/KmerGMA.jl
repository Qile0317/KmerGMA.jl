module KmerGMA
    using Plots, BioSequences, FASTX, Distances, #WebBlast, Statistics,
    ProgressMeter, #StringDistances, StatsBase, StatsPlots, DataFrames
    include("simpleExplore.jl")
    include("ExactMatch.jl")
    include("RefGen.jl")
    include("GMA.jl")
    include("API.jl")
    #include("GMA_Nless.jl")
end
