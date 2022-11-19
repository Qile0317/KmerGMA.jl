using KmerGMA
using Test, BioSequences, FASTX

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
    @test open(FASTX.FASTA.Reader,tf) do io
        recordCount(io) == 84
    end

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

    #testing GMA euclidean version
    kd = genKmers(1,withN=true)
    reference = open(FASTX.FASTA.Reader, tf)
    refKFV = kfv(genRef(1,reference,kd),kd) #generation of kmer frequency dict and converting
    close(reference)

    @test open(FASTX.FASTA.Reader,tf) do io
        eucGMA(k = 1, record = first(io), refVec = refKFV,
        windowsize = 289, kmerDict = kd, thr = 200.0, buff = 25,
        rv = fill(0.0,5)) == [71.2219387755102, 84.10289115646258,
        67.84098639455782, 67.84098639455782, 61.864795918367335,
        84.10289115646258, 95.00765306122447, 113.88860544217685]
    end

    #trying the testing version of the gma for a record
    kd = genKmers(6,withN = true)
    reference = open(FASTX.FASTA.Reader, tf)
    refKFV = kfv(genRef(6,reference,kd),kd) #generation of kmer frequency dict and converting
    close(reference)

    #def variables
    inp = FASTA.Record[]
    target = open(FASTX.FASTA.Reader, gf)
    goal = first(target)
    close(target)

    test_gma(k=6,record=goal,
    refVec = refKFV, windowsize = 289,
    kmerDict = kd,
    path = inp,
    thr = 250.0,
    buff = 20,
    rv= fill(0.0,5^6),
    thrbuff = " test ")

    """
    BenchmarkTools.Trial: 668 samples with 1 evaluation.
    Range (min … max):  6.143 ms … 12.836 ms  ┊ GC (min … max): 0.00% … 33.89%
    Time  (median):     6.983 ms              ┊ GC (median):    0.00%
    Time  (mean ± σ):   7.479 ms ±  1.620 ms  ┊ GC (mean ± σ):  8.96% ± 13.45%

       █▅     █▄▁
      ███▅▄▃▄████▄▃▃▂▂▃▂▂▂▁▁▁▁▁▁▁▁▁▁▁▂▁▁▁▁▁▁▁▂▃▄▄▃▂▂▂▃▅▄▃▃▂▂▂▁▂▂ ▃
      6.14 ms        Histogram: frequency by time        11.9 ms <

     Memory estimate: 7.67 MiB, allocs estimate: 164551.
    """

    @test inp == [FASTA.Record("AM773548.1 | SED = 98.17 | Pos = 6852:7141 test ",
    dna"GGTCCGTCAGG"), FASTA.Record("AM773548.1 | SED = 130.7 | Pos = 33845:34134 test ",
    dna"CAATGCCATGG"), FASTA.Record("AM773548.1 | SED = 249.1 | Pos = 33953:34242 test ",
    dna"CCATGGGCTGG")]

    @testset "gma_return_mode" begin
        #reference, genome is already deinfed
        res = FASTA.Record[FASTA.Record("test", dna"tt")]

        gma(k=6,record=goal, refVec=refKFV,windowsize =289,
        kmerDict = kd, thr = 200.0, buff = 50, rv=fill(0.0,5^6),
        thrbuff = "test", mode = "return", resultVec = res)

        """
        BenchmarkTools.Trial: 642 samples with 1 evaluation.
        Range (min … max):  6.190 ms … 17.251 ms  ┊ GC (min … max): 0.00% … 34.24%
        Time  (median):     7.066 ms              ┊ GC (median):    0.00%
        Time  (mean ± σ):   7.785 ms ±  1.866 ms  ┊ GC (mean ± σ):  9.10% ± 13.52%

        ▂▇▂  ▁█▄▂
        ███▅▃████▄▄▃▃▃▃▃▃▃▂▂▁▂▂▂▁▂▂▂▁▁▁▁▂▃▄▃▃▃▃▃▄▄▂▂▂▁▂▂▂▃▂▁▂▂▂▁▂▂ ▃
        6.19 ms        Histogram: frequency by time        13.5 ms <

        Memory estimate: 7.68 MiB, allocs estimate: 164519.
        """

        @test res[1] == FASTA.Record("test", dna"tt") #see if the first record stayed
        @test getSeq(res[2])[1:30] == dna"GTCTGGGGGAGGCTTGGTGCAGCCTGGGGG"
        @test getSeq(res[3])[1:30] == dna"CAGGCTCAGGTGCAGCTGGTGGAGTCTGGG" #not sure abt if the seq buffer worked..
        @test FASTA.description(res[3]) == "AM773548.1 | SED = 130.7 | Pos = 33845:34134test"
    end
end

@testset "API.jl" begin
    #try a version with a lowered SED
    @test testFindGenes(genome = GF, ref = tf, thr = 100.0) == [FASTX.FASTA.Record(
    "AM773548.1 | SED = 98.17 | Pos = 6852:7141 | thr = 100.0 | buffer = 50",
    dna"TCCTGACCAGG")]

    @test testFindGenes(genome = GF, ref = tf) == FASTX.FASTA.Record[FASTA.Record(
    "JQ684648.1 | SED = 130.8 | Pos = 8543:8832 | thr = 371.0 | buffer = 50", dna"CCGATTCACCA"),
    FASTA.Record("JQ684648.1 | SED = 110.6 | Pos = 20425:20714 | thr = 371.0 | buffer = 50",
    dna"GAATCCATGAA"), FASTA.Record("AM773729.1 | SED = 130.8 | Pos = 685:974 | thr = 371.0 | buffer = 50",
    dna"GGCCGATTCAC"), FASTA.Record("AM773729.1 | SED = 110.6 | Pos = 12791:13080 | thr = 371.0 | buffer = 50",
    dna"GAATCCATGAA"), FASTA.Record("AM773548.1 | SED = 98.17 | Pos = 6852:7141 | thr = 371.0 | buffer = 50",
    dna"GACTCCGTGAA"), FASTA.Record("AM773548.1 | SED = 368.9 | Pos = 23826:24115 | thr = 371.0 | buffer = 50",
    dna"TGCCTGGTGGC"), FASTA.Record("AM773548.1 | SED = 298.9 | Pos = 23931:24220 | thr = 371.0 | buffer = 50",
    dna"TGATGGCAGCA"), FASTA.Record("AM773548.1 | SED = 130.7 | Pos = 33845:34134 | thr = 371.0 | buffer = 50",
    dna"ATTCACCATCT")] #takes about 125 miliseconds

    #@test findGenes(genome = GF, ref = tf) == [] <- why doesnt it work...
end

@testset "ExactMatch.jl" begin
    #exactMatch set - single sequence
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
