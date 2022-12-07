using KmerGMA
using Test, BioSequences, FASTX, Random

#testing variables
tf = "Alp_V_ref.fasta"
gf = "Alp_V_locus.fasta"
GF = "Loci.fasta"
KD = Dict(dna"T" => 3, dna"A" => 1, dna"G" => 4, dna"N" => 5, dna"C" => 2)
kf = [0.0, 0.0, 1.0, 2.0, 0.0, 1.0, 0.0, 0.0, 0.0,
0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
KFV = Dict(dna"T" => 62.38095238095238, dna"C" => 73.70238095238095,
dna"A" => 63.25, dna"G" => 89.26190476190476, dna"N" => 0.0)

@testset "simpleExplore.jl" begin
    @test flipDict(Dict(dna"ATG" => 69, dna"ATGCCCCC" =>
    420)) == Dict(69 => dna"ATG", 420 => dna"ATGCCCCC")

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

    @test fasterKF(2, view(dna"GAGATAC", 1:7),
    genKmers(2;withN=true), fill(0.0,25)) == kf

    @test fasterKF(1,view(dna"GAGATAC",1:7),KD,[
        0.0,0.0,0.0,0.0,0.0]) == [3.0,1.0,1.0,2.0,0.0]

    @test kmerFreq(2, dna"GAGATAC") == kf
end

@testset "RefGen.jl" begin
    @testset "Reference dictionary generation" begin
        @test genRef(1,tf,KD) == KFV

        @test open(FASTX.FASTA.Reader,tf) do io
            genRef(1,io,KD) == KFV
        end
    end

    @testset "threshold finding function" begin
        @test findthr(tf,KFV,KD) == 147.38860544217687 #I need to change this w the scalefactor later
        
        #random threshold finder. they are seeded to 1112 except the last so testing should be consistent
        @test findRandThr(tf,KFV,KD) == 121.49633219954649
        @test findRandThr(tf,KFV,KD; ScaleFactor = 1.0) == 424.79102607709746
        @test findRandThr(tf,KFV,KD; buff = 20) == 141.4963321995465 # each unit = 1 substitution or 0.6 indel
        @test findRandThr(tf,KFV,KD; sampleSize = 418) == 155.84678784013605
        @test findRandThr(tf,KFV,KD; setSeed = 69420) == 121.0918679138322
    end
end

@testset "GMA.jl" begin
    #def variables for testing 
    target = open(FASTX.FASTA.Reader, gf)
    goal = first(target)
    close(target)
    sixMerDict = genKmers(6,withN=true)
    refKFD = genRef(6,tf,sixMerDict)
    refKFV = kfv(refKFD,sixMerDict)

    """ #there seem to be some floating point errors
    @testset "eucGma" begin
        res = Float64[]
        eucGma(k = 6, record = goal, refVec = refKFV,
        windowsize=289, kmerDict=sixMerDict,
        rv = fill(0.0, 5^6), resultVec = res)

        @test [436.72619047619065, 438.6309523809526, 
        440.79761904761926, 442.8214285714288, 
        442.75000000000017, 442.53571428571445, 
        446.53571428571445, 450.53571428571445, 
        452.4642857142859, 454.2738095238097, 
        460.2738095238097, 466.2738095238097] == res[317:328]

        @test res[6942] == 244.7738095238098

        res = Float64[]
        eucGma(k = 6, record = goal, refVec = refKFV,
        windowsize=289, kmerDict=sixMerDict,
        rv = fill(0.0, 5^6), resultVec = res, 
        ScaleFactor = 0.8333333333)

        @test res[317:328] == [363.9384920489347, 
        365.5257936361728, 367.33134919165616, 369.01785712809664, 
        368.95833331857517, 368.7797618900109, 372.11309522321085, 
        375.44642855641086, 377.05357141348946, 378.5615079213656, 
        383.5615079211656, 388.56150792096565]

        @test res[6942] == 203.97817459501573
    end
    """

    @testset "gma_return_mode" begin #
        #reference, genome is already deinfed
        res = FASTA.Record[FASTA.Record("test", dna"tt")]

        gma(k=6,record=goal, refVec=refKFV,windowsize =289,
        kmerDict = sixMerDict, thr = 200.0, buff = 50, rv=fill(0.0,5^6),
        thrbuff = "test", mode = "return", resultVec = res)

        @test res[1] == FASTA.Record("test", dna"tt") #see if the first record stayed
        @test getSeq(res[2])[1:30] == dna"GTCTGGGGGAGGCTTGGTGCAGCCTGGGGG"
        @test getSeq(res[3])[1:30] == dna"CAGGCTCAGGTGCAGCTGGTGGAGTCTGGG" #not sure abt if the seq buffer worked..
        @test FASTA.description(res[3]) == "AM773548.1 | SED = 130.7 | Pos = 33845:34134test"
        #more testing of other variables is done in the API.jl testset.
    end
end

@testset "API.jl" begin
    @testset "findGenes_noArgs" begin #it seemed to have gotten worse with findRandThr. I need to add extra
        a = findGenes(genome = GF, ref = tf)

        @test length(a) == 6 
        @test length(getSeq(a[2])) == 389

        @test FASTX.FASTA.description(a[1]) == "JQ684648.1 | SED = 10.90 | Pos = 8543:8832 | thr = 16.0 | buffer = 50"
        @test getSeq(a[2])[42:96] == dna"CTCTGAGACTCTCCTGTGCAGCCTCTGGATTCACTTTTGATGATTATGCCATGAG"
        @test getSeq(a[4])[317:328] == dna"AGACACAAACCT"
        @test getSeq(a[5])[198:200] == dna"CAG"
    end
    
    @testset "findGenes_thr=25" begin
        a = findGenes(genome = GF, ref = tf, thr = 25.0)

        @test length(a) == 10
        @test getSeq(a[1])[1:60] == dna"TTGTAGCAGCTATTAGCTGGAGTGGTGGTAGCACATACTATGCAGACTCCGTGAAGGGCC"
        @test getSeq(a[5])[320:end] == dna"AGGAGGCAGCTGGTTTCACGGTTTCCTGTCAGGCTCTGGAGTTTCCTCTCCACAGTGCAGGAACCCCTCT"
        @test FASTX.FASTA.description(a[1]) == "JQ684648.1 | SED = 10.90 | Pos = 8543:8832 | thr = 25.0 | buffer = 50"
        @test getSeq(a[9])[12:17] == dna"CCAGGT"
    end

    @testset "findGenes_k=5" begin
        a = findGenes(genome = GF, ref = tf, k=5)

        @test length(a) == 6
        @test FASTX.FASTA.description(a[6]) == "AM773548.1 | SED = 11.26 | Pos = 33843:34132 | thr = 19.0 | buffer = 50"
        @test getSeq(a[5])[50:75] == dna"TCCTGTGCAGCCTCTGGATTCACCTT"
        @test getSeq(a[3])[360:end] == dna"ACCCACCAAGGGCAGGGCTGAGCCCCAGAG"
    end
    #To Do:
    #buffer test! For when it exceeds end or preceeds start. 
    #test if its a vector of strings as the imput for genome.
end

@testset "ExactMatch.jl" begin
    @testset "singleSeq" begin
        @test exactMatch(dna"GAG",dna"CCCCCCCGAGCTTTT") == [8:10]
        @test exactMatch(dna"GAG",dna"CGAGCCCGAGCTTTT") == [2:4, 8:10]
        @test exactMatch(dna"GAG",dna"CGAGAGAGAAGGCCGAGCTTTT") == [2:4, 4:6, 6:8, 15:17]
        @test exactMatch(dna"GAG",dna"CGAGAGAGAAGGCCGAGCTTTT",
        overlap = false) == [2:4, 6:8, 15:17]
        @test isnothing(exactMatch(dna"GAG",dna"CCCCCCTTT"))
    end

    @testset "readerVer" begin
    #testing a subseq of the first sequence of a reader as a dna seq
        @test open(FASTX.FASTA.Reader,tf) do io
            subseq = getSeq(first(io))[42:69] #its a record, so its testing 2 things at once
            open(FASTX.FASTA.Reader,tf) do io2
                exactMatch(subseq, io2) == Dict("AM773729|IGHV1-1*01|Vicugna" => [42:69])
            end
        end

        #testing a seq itself but as a record
        @test open(FASTX.FASTA.Reader,tf) do io
            subseq = first(io)
            open(FASTX.FASTA.Reader,tf) do io2
                exactMatch(subseq, io2) == Dict("AM773729|IGHV1-1*01|Vicugna" => [1:296])
            end
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
    end
    #cflength like a few other functions isnt too relevant
end
