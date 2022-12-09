module KmerGMA
    using BioSequences, FASTX, Distances, Random
     #WebBlast, Statistics, StringDistances, StatsBase, StatsPlots, DataFrames
    include("simpleExplore.jl")
    include("ExactMatch.jl")
    #include("kmerCount.jl")
    include("GMAutils.jl")
    include("RefGen.jl")
    include("GMA.jl")
    include("GMA_Nless.jl")
    include("API.jl")
end