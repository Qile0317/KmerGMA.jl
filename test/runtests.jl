using KmerGMA
using Test, BioSequences, FASTX

#testing variables
tf = "Alp_V_ref.fasta"
gf = "Alp_V_locus.fasta"
KD = Dict(dna"T" => 3, dna"A" => 1, dna"G" => 4, dna"N" => 5, dna"C" => 2)
kf = [0.0, 0.0, 1.0, 2.0, 0.0, 1.0, 0.0, 0.0, 0.0,
0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
KFV = Dict(dna"T" => 62.38095238095238, dna"C" => 73.70238095238095,
dna"A" => 63.25, dna"G" => 89.26190476190476, dna"N" => 0.0)

@testset "simpleExplore.jl" begin
    @test open(FASTX.FASTA.Reader,tf) do io
        recordCount(io) == 84
    end

    @test flipDict(Dict(dna"ATG" => 69, dna"ATGCCCCC" => 420)) == Dict(69 => dna"ATG", 420 => dna"ATGCCCCC")

    @test kfv(Dict(dna"A" => 69.0, dna"G" => 420.0,
    dna"C" => 0.0,dna"N" => 1.0, dna"T" => 42.0),
    genKmers(1; withN=true)) == [69.0, 0.0, 42.0, 420.0, 1.0]

    @test percentN(dna"ACGN") == 0.25

    @test avgRecLen(tf) == 289

    @test open(FASTX.FASTA.Reader,tf) do io
        readerNTs(io) == 24242
    end
end

@testset "GMAutils.jl" begin
    @test open(FASTX.FASTA.Reader,tf) do io
        getSeq(first(io))[10:20] == dna"CTGGTGCAGCC"
    end

    @test genKmers(1; withN=true) == KD

    @test fasterKF(2, dna"GAGATAC",
    genKmers(2;withN=true), fill(0.0,25)) == kf

    @test kmerFreq(2,dna"GAGATAC") == kf
end

@testset "RefGen.jl" begin
    @test genRef(1,tf,KD) == KFV

    @test open(FASTX.FASTA.Reader,tf) do io
        genRef(1,io,KD) == KFV
    end

    @test open(FASTX.FASTA.Reader,tf) do io
        findthr(io,KFV,KD) == 147.38860544217687
    end

    @test findthr(tf,KFV,KD) == 147.38860544217687
end

@testset "GMA.jl" begin
    @test fasterKF(1,dna"GAGATAC",KD,[0.0,0.0,0.0,0.0,0.0]) == [3.0,1.0,1.0,2.0,0.0]

    @test open(FASTX.FASTA.Reader,tf) do io
        queryMatch(1,first(io),KFV,KD,289) == [71.2219387755102,
        84.10289115646258, 84.10289115646258, 84.10289115646258, 80.1267006802721,
        104.36479591836734, 117.26955782312923]
    end #this is outdated and probably doesnt work.

    #@test open(FASTX.FASTA.Reader,tf) do reference
    #    open(FASTX.FASTA.Reader,gf) do target
    #        testFindGenes(genome=target,ref=reference) == "placeholder"
    #    end
    #end

    #problem: now that im using @views, i may have to revamp everything to the type of seqView...
end

@testset "API.jl" begin
    reference = open(FASTX.FASTA.Reader, tf)
    target = open(FASTX.FASTA.Reader, gf)
    #def variables
    kd = genKmers(6,withN=true)
    k = 6
    refKFD = genRef(k,reference,kd) #generation of kmer frequency dict
    refKFV = kfv(refKFD,kd)
    inp = FASTA.Record[]

    test_gma(k=k,record=first(target),
    refVec = refKFV, windowsize = 289,
    kmerDict = kd,
    path = inp,
    thr = 250.0,
    buff = 20,
    rv= fill(0.0,5^6),
    thrbuff = " test ")

    @test inp == [FASTA.Record("AM773548.1 | SED = 98.17 | Pos = 6852:7141 test ",
    dna"GGTCCGTCAGG"), FASTA.Record("AM773548.1 | SED = 130.7 | Pos = 33845:34134 test ",
    dna"CAATGCCATGG"), FASTA.Record("AM773548.1 | SED = 249.1 | Pos = 33953:34242 test ",
    dna"CCATGGGCTGG")]

    close(reference)
    close(target)

    #testing the test api, but the reading in of genomes is bugged atm. idk why reading in locus yields no length
    #@test open(FASTX.FASTA.Reader,tf) do reference
    #    open(FASTX.FASTA.Reader,gf) do target
    #        testFindGenes(genome = target, ref = reference)
    #    end
    #end
end

@testset "ExactMatch.jl" begin
    #exactMatch set - single sequence
    @test exactMatch(dna"GAG",dna"CCCCCCCGAGCTTTT") == [8:10]
    @test exactMatch(dna"GAG",dna"CGAGCCCGAGCTTTT") == [2:4, 8:10]
    @test exactMatch(dna"GAG",dna"CGAGAGAGAAGGCCGAGCTTTT") == [2:4, 4:6, 6:8, 15:17]
    @test exactMatch(dna"GAG",dna"CGAGAGAGAAGGCCGAGCTTTT",
    overlap = false) == [2:4, 6:8, 15:17]
    @test exactMatch(dna"GAG",dna"CCCCCCTTT") == nothing
    @test open(FASTX.FASTA.Reader,tf) do io
        exactMatch(FASTA.sequence(first(io))))[42:69],
        io) == Dict("AM773729|IGHV1-1*01|Vicugna" => [42:69])
    end

    #exactMatch - reader
    @test open(FASTX.FASTA.Reader,tf) do io
        exactMatch(dna"AAAAAAAAA", io) == "no match"
    end
    @test open(FASTX.FASTA.Reader,tf) do io
        exactMatch(dna"AAATT", io) ==
        Dict("AM773729|IGHV1-1*01|Vicugna" => [174:178],
        "AM939700|IGHV1S5*01|Vicugna" => [174:178])
    end
    @test open(FASTX.FASTA.Reader,tf) do io
        exactMatch(FASTA.sequence(first(io)),
        io) == Dict("AM773729|IGHV1-1*01|Vicugna" => [1:296])

    #cflength like a few other functions depend on FASTA.seqlen() which doesnt work in testing for some reason...
end
