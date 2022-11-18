module KmerGMA
    using BioSequences, FASTX, Distances, TranscodingStreams
     #WebBlast, Statistics, StringDistances, StatsBase, StatsPlots, DataFrames
    include("simpleExplore.jl")
    include("ExactMatch.jl")
    include("GMAutils.jl")
    include("RefGen.jl")
    include("GMA.jl")
    include("API.jl")
end
