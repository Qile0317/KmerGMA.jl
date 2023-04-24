using Test, BioSequences, FASTX, Random, BioAlignments, Distances, StaticArrays

# for devs: if testing the script on your own machine, set the following variable to true. Otherwise must set to false when pushing
test_locally = false

# testing fasta files
tf = "fasta_files/Alp_V_ref.fasta"
test_mini_genome = "Alp_V_locus.fasta"
test_genome = "Loci.fasta"
test_8_seqs = "fasta_files/8_ident_Alp_V_loci.fasta"

# load the package
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

    # strobemer stuff
    include("../src/StrobemerGMA/Strobemers.jl")
    include("../src/StrobemerGMA/StrobeRefGen.jl")
    include("../src/StrobemerGMA/MonteCarloBenchmark.jl")
    include("../src/StrobemerGMA/StrobeGenomeMiner.jl")

    # multithreading
    include("../src/MultiThread/GenomeMiner.jl")

    tf = "test/" * tf
    test_mini_genome = "test/" * test_mini_genome
    test_genome = "test/" * test_genome
    test_8_seqs = "test/" * test_8_seqs
else
    using KmerGMA
end

# testing variables
const test_seq = dna"ATGCATGC"
const test_consensus_seq = dna"CAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTGCAGCCTGGGGGGTCTCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTATGCCATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTCGAGTGGGTCTCAGCTATTAATAGTGGTGGTGGTAGCACATACTATGCAGACTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAAACCTGAGGGCACGGCCGTGTATTACTGTGGTAAAGAAGA"

const test_KFV = [0.0, 0.0, 0.0, 2.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0]

# mean function surprising doesnt exist in base library
function mean(vec::Vector{Float64})
    len, curr_sum = 0, 0
    for num in vec
        curr_sum += num
        len += 1
    end
    return curr_sum/len
end

# run tests 
println("%%%%%%%%%% KmerGMA %%%%%%%%%%")
include("test_folder/test-KmerGMA.jl")

println("%%%%%%%%%% StrobemerGMA %%%%%%%%%%")
include("test_folder/test-StrobemerGMA.jl")