using KmerGMA
using Test, BioSequences, FASTX, Random, BioAlignments

#testing variables
const tf = "Alp_V_ref.fasta"
const test_mini_genome = "Alp_V_locus.fasta"
const test_genome = "Loci.fasta"

# testing sequences 
const test_seq = dna"ATGCATGC"
const test_consensus_seq = dna"CAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTGCAGCCTGGGGGGTCTCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTATGCCATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTCGAGTGGGTCTCAGCTATTAATAGTGGTGGTGGTAGCACATACTATGCAGACTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAAACCTGAGGGCACGGCCGTGTATTACTGTGGTAAAGAAGA"

const test_KFV = [0.0, 0.0, 0.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0]

# mean function surprising doesnt exist in base library
mean(vec::Vector{Float64}) = sum(vec)/length(vec)

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

    @testset "kmer bit utils" begin
        @test as_UInt(test_seq) == unsigned(14649)
        @test as_kmer(unsigned(14649), 8) == test_seq
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
        record_vec = nothing
    end

    @testset "cluster_ref_API" begin
        @test get_cluster_index(5, [1,2,6,10]) == 3
        @test get_cluster_index(12, [1,2,6,10]) == 5
        @test get_cluster_index(0, [1,2,6,10]) == 1
        
        @testset "include_avg = false" begin
            a = cluster_ref_API(tf, 1, cutoffs = [7,12,20,25], include_avg = false)
            @test a[1] == [[62.785714285714285, 72.78571428571429, 89.78571428571429, 62.642857142857146], [63.13333333333333, 71.33333333333333, 90.53333333333333, 62.6], [63.5, 70.71428571428571, 90.78571428571429, 64.07142857142857], [62.54545454545455, 68.72727272727273, 91.36363636363636, 64.54545454545455], [63.666666666666664, 78.53333333333333, 86.9, 60.56666666666667]]
            @test a[2] == [288,288,289,287,290]
            @test length(a[3]) == 5
            @test a[3][1][1:4] == dna"CAGG" # doesen't test all sequences lol but they should be fine if all prev tests passed
            @test a[4] == [false,false,false,false,false]
            a = nothing # save mem
        end

        @testset "include_avg = true (default)" begin
            a = cluster_ref_API(tf, 1, cutoffs = [7,12,20,25])
            @test a[1] == [[62.785714285714285, 72.78571428571429, 89.78571428571429, 62.642857142857146], [63.13333333333333, 71.33333333333333, 90.53333333333333, 62.6], [63.5, 70.71428571428571, 90.78571428571429, 64.07142857142857], [62.54545454545455, 68.72727272727273, 91.36363636363636, 64.54545454545455], [63.666666666666664, 78.53333333333333, 86.9, 60.56666666666667], [63.25, 73.70238095238095, 89.26190476190476, 62.38095238095238]]
            @test a[2] == [288, 288, 289, 287, 290, 289]
            @test length(a[3]) == 6
            @test a[3][1][1:4] == dna"CAGG" # doesen't test all sequences lol but they should be fine if all prev tests passed
            @test a[4] == [false,false,false,false,false, false]
            a = nothing # save mem
        end
    end

    @testset "eliminate_null_params" begin
        test_KFVs = Vector{Float64}[test_KFV,test_KFV.+0.3]
        test_ws, test_cons_vec = Int[8,9], Seq[test_seq, test_seq*dna"Y"]
        inv_vec = [false,true]

        test_KFVs, test_ws, test_cons_vec = eliminate_null_params(test_KFVs, test_ws, test_cons_vec, inv_vec)

        @test test_KFVs == [test_KFV]
        @test test_ws == [8]
        @test test_cons_vec == [test_seq]

        RVs, windowsizes, consensus_refseqs, invalids = cluster_ref_API(tf, 6; cutoffs = [7,12,20,25])
        RVs, windowsizes, consensus_refseqs = eliminate_null_params(RVs, windowsizes, consensus_refseqs, invalids)
        @test windowsizes == [288,288,288,289,290,289]
        @test length(RVs) == length(consensus_refseqs) == 6
    end
end

@testset "DistanceTesting.jl" begin
    RV, ws, cons_seq = gen_ref_ws_cons(tf, 6)
    @test Int(round(estimate_optimal_threshold(RV,299; buffer = 12))) == 27
    
    rvs,ws,cons,inv = cluster_ref_API(tf, 6; cutoffs = [7,12,20,25], include_avg = false)
    rounded_res = [Int(round(num)) for num in estimate_optimal_threshold(rvs,ws; buffer = 8)]
    @test rounded_res == [38, 33, 41, 37, 29]
end

@testset "Alignment.jl" begin
    score_model = AffineGapScoreModel(EDNAFULL, gap_open=-5, gap_extend=-1)
    aln_obj = pairalign(SemiGlobalAlignment(),dna"ATGCATGC",dna"GGGGGATGCATGCAAAAA",score_model)
    @test cigar_to_UnitRange(aln_obj) == 6:13

    aln_obj = pairalign(SemiGlobalAlignment(),dna"ATGCATGC",dna"GGGGGATGCTTATGCAAAAA",score_model)
    @test cigar_to_UnitRange(aln_obj) == 6:15
end

@testset "RSS.jl" begin # UNFINISHED
    rss_align = Align_RSS(test_seq*HumanRSSV*test_seq, HumanRSSV)
    @test cigar(rss_align.aln.a.aln) == "8D28=8D"

    @test RSS_dist(HumanRSSV, HumanRSSV) == 0
    @test RSS_dist(HumanRSSV[1:end-1]*dna"T", HumanRSSV) == 1

    @test is_RSS(rss_align) == true
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

        @test round(mean(dist_vec)) == 46

        @test length(res) == 3
        @test getSeq(res[1]) == getSeq(res[2])
        @test FASTA.description(res[3]) == "AM773548.1 | dist = 8.18 | MatchPos = 6801:7189 | GenomePos = 444023"
    end
end

@testset "OmnGenomeMiner.jl" begin #unfinished
    @testset "buff = 200" begin
        rvs,ws,cons,inv = cluster_ref_API(tf, 6; cutoffs = [7,12,20,25], include_avg = false)
        test_res = FASTA.Record[]

        # thresholds are purposefully lower
        Omn_KmerGMA!(genome_path = test_mini_genome, refVecs = rvs, windowsizes = ws, consensus_seqs = cons, resultVec = test_res,
                    buff = 200, thr_vec = [37,33,38,34,28,27])
        @test length(test_res) == 5
        @test FASTA.description(test_res[5]) == "AM773548.1 | Dist = 33.1 | KFV = 2 | MatchPos = 33845:34132 | GenomePos = 0 | Len = 287"
        @test FASTA.description(test_res[4]) == "AM773548.1 | Dist = 37.11 | KFV = 1 | MatchPos = 33845:34132 | GenomePos = 0 | Len = 287"
        @test FASTA.description(test_res[3]) == "AM773548.1 | Dist = 38.0 | KFV = 3 | MatchPos = 33845:34132 | GenomePos = 0 | Len = 287"
        @test FASTA.description(test_res[2]) == "AM773548.1 | Dist = 34.12 | KFV = 4 | MatchPos = 23907:24179 | GenomePos = 0 | Len = 272"
        @test FASTA.description(test_res[1]) == "AM773548.1 | Dist = 38.0 | KFV = 3 | MatchPos = 6852:7139 | GenomePos = 0 | Len = 287"
    end
end

@testset "API.jl" begin
    a = findGenes(genome_path = test_mini_genome, ref_path = tf)[1]
    @test length(a) == 3
    @test FASTA.description(a[1]) == "AM773548.1 | dist = 8.18 | MatchPos = 6852:7153 | GenomePos = 0"
    @test FASTA.description(a[2]) == "AM773548.1 | dist = 24.91 | MatchPos = 23907:24249 | GenomePos = 0"
    @test FASTA.description(a[3]) == "AM773548.1 | dist = 10.9 | MatchPos = 33845:34144 | GenomePos = 0"

    # more comprehensive testing is needed of other params but my personal testing shows that it works nicely
    # need to test findGenes_cluster_mode
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
    
    @testset "fasta_id_to_cumulative_len_dict" begin
        @test fasta_id_to_cumulative_len_dict(test_genome) == Dict{String, Int64}(
            "JQ684648.1 Lama glama clone V03 IgH locus genomic sequence" => 0,
            "JQ684647.1 Lama glama clone F07 IgH locus genomic sequence" => 121478,
            "AM773548.1 Lama pacos germline IgHV region, Vh3-S1, Vh2-S1 and vhh3-S1 genes" => 444023,
            "AM773729.1 Lama pacos germline IgH locus: proximal IgHV region genes, complete IgHD region genes, complete IgHJ region genes and complete IgHC region genes" => 221227
            )
    end
end

# work in progress paired kmers for upcoming features
@testset "PairedKmers.jl" begin
    @test (unsigned(228), unsigned(228)) == initialize_kmers(test_seq, 6)
    @test dna"AT" == as_kmer(initialize_kmers(test_seq, 3)[1], 2)

    @test Int(as_index(as_UInt(dna"ATGC"), as_UInt(dna"ATGC"), 4)) == 14650

    @test kmer_pair_count(test_seq, 1) == fill(4.0, 16)
    @test sum(kmer_pair_count(test_seq, 2)) == 49.0
    @test round(mean(kmer_pair_count(test_seq, 2)), digits = 5) == 0.19141

    @testset "kmer_pair_count!" begin
        bins = zeros(16)
        kmer_pair_count!(test_seq,1,8,bins,unsigned(3))
        @test bins == fill(4.0, 16)

        bins = zeros(256)
        kmer_pair_count!(test_seq,2,8,bins,unsigned(63))
        @test sum(bins) == 49.0
        @test round(mean(bins), digits = 5) == 0.19141
    end
end