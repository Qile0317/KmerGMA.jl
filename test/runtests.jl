using Test, BioSequences, FASTX, Random, BioAlignments, Distances, StaticArrays

# for devs: if testing the script on your own machine, set the following variable to true. Otherwise must set to false when pushing
test_locally = false

#setting testing variables
if test_locally
    include("../src/Consts.jl")
    include("../src/Kmers.jl")
    include("../src/Consensus.jl")
    include("../src/DistanceTesting.jl")
    include("../src/ReferenceGeneration.jl")
    include("../src/Alignment.jl")
    include("../src/GenomeMiner.jl")
    include("../src/OmnGenomeMiner.jl")
    include("../src/API.jl")
    include("../src/ExactMatch.jl")

    # in progress
    include("../src/RSS.jl")
    include("../src/PairedKmers.jl")
    include("../src/StrobemerGMA/Strobemers.jl")
    include("../src/StrobemerGMA/MonteCarloBenchmark.jl")
    include("../src/StrobemerGMA/StrobeGenomeMiner.jl")

    tf = "test/Alp_V_ref.fasta"
    test_mini_genome = "test/Alp_V_locus.fasta"
    test_genome = "test/Loci.fasta"
else
    using KmerGMA
    tf = "Alp_V_ref.fasta"
    test_mini_genome = "Alp_V_locus.fasta"
    test_genome = "Loci.fasta"
end

# testing sequences 
const test_seq = dna"ATGCATGC"
const test_consensus_seq = dna"CAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTGCAGCCTGGGGGGTCTCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTATGCCATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTCGAGTGGGTCTCAGCTATTAATAGTGGTGGTGGTAGCACATACTATGCAGACTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAAACCTGAGGGCACGGCCGTGTATTACTGTGGTAAAGAAGA"

const test_KFV = [0.0, 0.0, 0.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0]

# mean function surprising doesnt exist in base library
mean(vec::Vector{Float64}) = sum(vec)/length(vec)

# test 
include("test_folder/test-KmerGMA.jl")
include("test_folder/test-StrobemerGMA.jl")