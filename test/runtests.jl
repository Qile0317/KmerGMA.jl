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
    include("../src/StrobemerGMA/StrobeRefGen.jl")
    include("../src/StrobemerGMA/MonteCarloBenchmark.jl")
    include("../src/StrobemerGMA/StrobeGenomeMiner.jl")

    # in progress multithreading
    include("../src/Multithreaded/GenomeMiner.jl")

    tf = "test/Alp_V_ref.fasta"
    test_mini_genome = "test/Alp_V_locus.fasta"
    test_genome = "test/Loci.fasta"
else
    using KmerGMA
    tf = "Alp_V_ref.fasta"
    test_mini_genome = "Alp_V_locus.fasta"
    test_genome = "Loci.fasta"
end

const test_seq = dna"ATGCATGC"
const test_consensus_seq = dna"CAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTGCAGCCTGGGGGGTCTCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTATGCCATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTCGAGTGGGTCTCAGCTATTAATAGTGGTGGTGGTAGCACATACTATGCAGACTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAAACCTGAGGGCACGGCCGTGTATTACTGTGGTAAAGAAGA"

const test_KFV = [0.0, 0.0, 0.0, 2.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0]

# mean function surprising doesnt exist in base library
function mean(vec::Vector{Float64})
    len::Int, curr_sum = 0, 0
    for num in vec
        curr_sum += num
        len += 1
    end
    return curr_sum/len
end

# test 
println("%%%%%%%%%% KmerGMA %%%%%%%%%%")
include("test_folder/test-KmerGMA.jl")

println("%%%%%%%%%% StrobemerGMA %%%%%%%%%%")
include("test_folder/test-StrobemerGMA.jl")