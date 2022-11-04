using KmerGMA
using Test

#testing variables
tf = "test/Alp_V_ref.fasta"
KD = Dict(dna"T" => 3, dna"A" => 1, dna"G" => 4, dna"N" => 5, dna"C" => 2)
kf = [0.0, 0.0, 1.0, 2.0, 0.0, 1.0, 0.0, 0.0, 0.0,
0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
KFV = Dict(dna"T" => 62.38095238095238, dna"C" => 73.70238095238095,
dna"A" => 63.25, dna"G" => 89.26190476190476, dna"N" => 0.0)

@testset "simpleExplore.jl" begin
    @test genKmers(1; withN=true) == KD

    @test kmerFreq(2,dna"GAGATAC") == kf

    @test recordCount(open(FASTX.FASTA.Reader,tf)) == 84

    @test avgRecLen(tf) == 289

    @test flipDict(Dict(dna"ATG" => 69, dna"ATGCCCCC" => 420)) == Dict(69 => dna"ATG", 420 => dna"ATGCCCCC")

    @test kfv(Dict(dna"A" => 69.0, dna"G" => 420.0,
    dna"C" => 0.0,dna"N" => 1.0, dna"T" => 42.0),
    genKmers(1; withN=true)) == [69.0, 0.0, 42.0, 420.0, 1.0]

    @test readerNTs(open(FASTX.FASTA.Reader,tf)) == 24242

    @test percentN(dna"ACGN") == 0.25
end

@testset "RefGen.jl" begin
    @test genRef(1,tf,KD) == KFV
    @test findthr(tf,KFV,KD) == 147.38860544217687
end

@testset "GMA.jl" begin
    @test fasterKF(1,dna"GAGATAC",KD,[0.0,0.0,0.0,0.0,0.0]) == [3.0,1.0,1.0,2.0,0.0]
    @test queryMatch(1,first(open(FASTX.FASTA.Reader,tf)),KFV,KD,289) == [71.2219387755102,
    84.10289115646258, 84.10289115646258, 84.10289115646258, 80.1267006802721,
    104.36479591836734, 117.26955782312923]
    #problem: the main function for writing cant be tested lol
    #problem: Blast shouldne be tested otherwise it'll destroy NCBI and get you banned
end
