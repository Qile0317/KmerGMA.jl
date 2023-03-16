using KmerGMA
using Test, BioSequences, FASTX, Random, BioAlignments

#testing variables
const tf = "Alp_V_ref.fasta"
const test_mini_genome = "Alp_V_locus.fasta"
const test_genome = "Loci.fasta"

# testing sequences 
const test_seq = dna"ATGCATGC"
const test_consensus_seq = dna"CAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTGCAGCCTGGGGGGTCTCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTATGCCATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTCGAGTGGGTCTCAGCTATTAATAGTGGTGGTGGTAGCACATACTATGCAGACTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAAACCTGAGGGCACGGCCGTGTATTACTGTGGTAAAGAAGA"

test_KFV = [0.0, 0.0, 0.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0]

@testset "Kmers.jl" begin
    @testset "kmer_count!" begin
        test_bins = zeros(4)
        kmer_count!(str = test_seq, k = 1, bins = test_bins, mask = unsigned(3))
        @test test_bins == [2,2,2,2]

        test_bins = zeros(16)
        kmer_count!(str = view(test_seq, 1:8), k = 2, bins = test_bins, mask = unsigned(15))
        @test test_bins == test_KFV
    end

    @testset "kmer_count" begin
        @test kmer_count(test_seq, 1) == [2,2,2,2]
        @test kmer_count(view(test_seq, 1:8), 2) == test_KFV
    end

    @testset "kmer_dist" begin
        @test kmer_dist(test_seq^25*dna"A"*test_seq^25, test_seq^25*dna"G"*test_seq^25, 2) == 1.0
        @test kmer_dist(test_seq^25*dna"AA"*test_seq^25, test_seq^25*dna"GT"*test_seq^25, 2) == 2.0
    end
end

@testset "Consensus.jl" begin
    @test Profile(2).vecs == [[0,0],[0,0],[0,0],[0,0]]
    @test Profile(3)[DNA_A] == [0,0,0]

    a = Profile(8)
    add_consensus!(a, test_seq)
    @test a.vecs == [[1, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 1], [0, 0, 1, 0, 0, 0, 1, 0], [0, 1, 0, 0, 0, 1, 0, 0]]
    @test a.len == 8

    lengthen!(a, 9)
    @test a.vecs == a.vecs == [[1, 0, 0, 0, 1, 0, 0, 0,0], [0, 0, 0, 1, 0, 0, 0, 1,0], [0, 0, 1, 0, 0, 0, 1, 0,0], [0, 1, 0, 0, 0, 1, 0, 0,0]]
    @test a.len == 9

    add_consensus!(a, test_seq[1:7]*dna"G"); add_consensus!(a, test_seq[1:7]*dna"G")
    @test a.vecs == [[3, 0, 0, 0, 3, 0, 0, 0, 0], [0, 0, 0, 3, 0, 0, 0, 1, 0], [0, 0, 3, 0, 0, 0, 3, 2, 0], [0, 3, 0, 0, 0, 3, 0, 0, 0]]
    @test a.len == 9

    @test consensus_seq(a)[1:8] == test_seq[1:7]*dna"G"
end

@testset "DistanceTesting.jl" begin
    @test estimate_optimal_threshold(RV,299) == 30.52158730158728
end

@testset "ReferenceGeneration.jl" begin
    @testset "gen_ref_ws_cons" begin
        @test gen_ref_ws_cons(tf, 1) == ([63.25, 73.70238095238095, 89.26190476190476, 62.38095238095238], 289, test_consensus_seq)
        @test gen_ref_ws_cons(tf, 1; get_maxlen = true) == ([63.25, 73.70238095238095, 89.26190476190476, 62.38095238095238], 289, test_consensus_seq, 299)
        @test gen_ref_ws_cons(tf,2)[1] == [11.178571428571429, 15.964285714285714, 
        24.154761904761905, 11.88095238095238, 22.76190476190476, 
        17.904761904761905, 8.154761904761905, 24.88095238095238, 
        18.607142857142858, 22.202380952380953, 30.369047619047617, 
        18.07142857142857, 10.702380952380953, 17.047619047619047, 
        26.166666666666664, 7.5476190476190474]

        @test gen_ref_ws_cons(tf,6)[1][5:10] == [0.011904761904761904, 
        0.023809523809523808, 0.0, 0.0, 0.023809523809523808, 0.0]

        record_vec, reader = FASTA.Record[], open(FASTA.Reader, tf)
        for record in reader; push!(record_vec, record) end; close(reader)

        @test gen_ref_ws_cons(record_vec, 1) ==  ([63.25, 73.70238095238095, 89.26190476190476, 62.38095238095238], 289, test_consensus_seq)
        @test gen_ref_ws_cons(record_vec, 1; get_maxlen = true) == ([63.25, 73.70238095238095, 89.26190476190476, 62.38095238095238], 289, test_consensus_seq, 299)
    end

    @testset "cluster_ref_API" begin
        @test get_cluster_index(5, [1,2,6,10]) == 3
        @test get_cluster_index(12, [1,2,6,10]) == 5
        @test get_cluster_index(0, [1,2,6,10]) == 1
        
        a = cluster_ref_API(tf, 1)
        @test a[1] == [[62.935483870967744, 71.90322580645162, 90.19354838709677, 62.774193548387096], [63.666666666666664, 70.83333333333333, 90.83333333333333, 63.916666666666664], [63.0, 69.07692307692308, 90.46153846153847, 64.76923076923077], [63.535714285714285, 79.07142857142857, 87.0, 60.17857142857143]]
        @test a[2] == [288,289,287,290]
        @test a[3][1][1:4] == dna"CAGG" # doesen't test all sequences lol but they should be fine if all prev tests passed
        @test a[4] == [false,false,false,false]
    end
end

@testset "Alignment.jl" begin
    score_model = AffineGapScoreModel(EDNAFULL, gap_open=-5, gap_extend=-1)
    aln_obj = pairalign(SemiGlobalAlignment(),dna"ATGCATGC",dna"GGGGGATGCATGCAAAAA",score_model)
    @test cigar_to_UnitRange(aln_obj) == 6:13

    aln_obj = pairalign(SemiGlobalAlignment(),dna"ATGCATGC",dna"GGGGGATGCTTATGCAAAAA",score_model)
    @test cigar_to_UnitRange(aln_obj) == 6:15
end

@testset "GenomeMiner.jl" begin
    RV, ws, cons_seq = gen_ref_ws_cons(tf, 6)

    @testset "do_align = false" begin
        consts = init_InputConsts(genome_path = test_genome,
        refVec = RV, consensus_refseq = cons_seq, windowsize = ws,
        thr = 30, do_align = false)

        res = FASTA.Record[FASTA.Record("test", dna"tt")]
        ac_gma_testing!(inp = consts, resultVec = res)

        @test length(res) == 8
        @test FASTA.description(res[3]) == "JQ684648.1 | dist = 30.0 | MatchPos = 20131:20519 | GenomePos = 0"
        @test FASTA.description(res[end-2]) == "AM773548.1 | dist = 8.18 | MatchPos = 6801:7189 | GenomePos = 444023"
    end

    @testset "do_align = true, get_hit_loci = true" begin
        res = FASTA.Record[]
        hit_vec = Int[]
        consts = init_InputConsts(genome_path = test_genome,
            refVec = RV, consensus_refseq = cons_seq, windowsize = ws,
            thr = 30, do_align = true,get_hit_loci = true)

        ac_gma_testing!(inp = consts, resultVec = res, hit_loci_vec = hit_vec)

        @test length(res) == 7
        @test hit_vec == [8543, 20175, 221912, 234016, 450875, 467930, 477868]
        @test FASTA.description(res[2]) == "JQ684648.1 | dist = 30.0 | MatchPos = 20175:20574 | GenomePos = 0"
        @test FASTA.description(res[end-2]) == "AM773548.1 | dist = 8.18 | MatchPos = 6852:7153 | GenomePos = 444023"
    end

    @testset "do_return_dists = true, buff = 0, thr = 10, do_align = false" begin
        res = FASTA.Record[]
        dist_vec = Float64[]
        consts = init_InputConsts(genome_path = test_genome,
            refVec = RV, consensus_refseq = cons_seq, windowsize = ws,
            thr = 10, do_align = false, do_return_dists = true)

        ac_gma_testing!(inp = consts, resultVec = res, dist_vec = dist_vec)

        @test length(dist_vec) == 484127
        @test mean(dist_vec) == 46.290337910562926

        @test length(res) == 3
        @test getSeq(res[1]) == getSeq(res[2])
        @test FASTA.description(res[3]) == "AM773548.1 | dist = 8.18 | MatchPos = 6801:7189 | GenomePos = 444023"
    end
end

@testset "API.jl" begin
    a = findGenes(genome_path = test_mini_genome, ref_path = tf)[1]
    @test length(a) == 3
    @test FASTA.description(a[1]) == "AM773548.1 | dist = 8.18 | MatchPos = 6852:7153 | GenomePos = 0"
    @test FASTA.description(a[2]) == "AM773548.1 | dist = 24.91 | MatchPos = 23907:24249 | GenomePos = 0"
    @test FASTA.description(a[3]) == "AM773548.1 | dist = 10.9 | MatchPos = 33845:34144 | GenomePos = 0"

    # more comprehensive testing is needed of other params but my personal testing shows that it works nicely
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