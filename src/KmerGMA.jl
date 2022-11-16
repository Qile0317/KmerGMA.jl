module KmerGMA
    using Plots, BioSequences, FASTX, Distances, ProgressMeter, StringViews #WebBlast, Statistics, StringDistances, StatsBase, StatsPlots, DataFrames
    include("simpleExplore.jl")
    include("ExactMatch.jl")
    include("GMAutils.jl")
    include("RefGen.jl")
    include("GMA.jl")
    include("API.jl")
    #include("BLAST.jl")
end
