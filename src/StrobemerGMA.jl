# new work in progress module
module StrobemerGMA
    using BioSequences, FASTX, Distances, Random, BioAlignments
    include("Strobemers.jl")
    include("StrobemerGMA/MonteCarloBenchmark.jl")
    include("StrobemerGMA/StrobeGenomeMiner.jl") # work in progress
end