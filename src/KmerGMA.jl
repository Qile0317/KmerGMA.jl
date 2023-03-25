module KmerGMA
    using BioSequences, FASTX, Distances, Random, BioAlignments, StaticArrays
    include("Consts.jl")
    include("Kmers.jl")
    include("Consensus.jl")
    include("DistanceTesting.jl")
    include("ReferenceGeneration.jl")
    include("Alignment.jl")
    include("GenomeMiner.jl")
    include("OmnGenomeMiner.jl")
    include("API.jl")
    include("ExactMatch.jl")

    # work in progress 
    include("RSS.jl")
    include("PairedKmers.jl")

    # strobemerGMA
    include("StrobemerGMA/Strobemers.jl")
    include("StrobemerGMA/MonteCarloBenchmark.jl")
    include("StrobemerGMA/StrobeGenomeMiner.jl")
end

"""
# new work in progress module
module StrobemerGMA
    using BioSequences, FASTX, Distances, Random, BioAlignments, StaticArrays
    include("Consts.jl")
    include("Kmers.jl")
    include("DistanceTesting.jl")

    include("StrobemerGMA/Strobemers.jl")
    include("StrobemerGMA/MonteCarloBenchmark.jl")
    include("StrobemerGMA/StrobeGenomeMiner.jl") # work in progress
end
"""