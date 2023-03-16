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

# immediate upcoming features for 0.5.0 or potentially 1.0.0:
# - improved threshold estimation via monte carlo testing. (package is currently optimized for k = 6)
# - Optional clustering of references, and the O(mn) KmerGMA