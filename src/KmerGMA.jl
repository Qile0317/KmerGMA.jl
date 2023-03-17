module KmerGMA
    using BioSequences, FASTX, Distances, Random, BioAlignments
    include("Consts.jl")
    include("Kmers.jl")
    include("Consensus.jl")
    include("DistanceTesting.jl")
    include("ReferenceGeneration.jl")
    include("Alignment.jl")
    include("GenomeMiner.jl")
    include("API.jl")
    include("ExactMatch.jl")
end