@testset "Strobemers.jl" begin
    @test randstrobe_score(dna"ATGC", dna"GTGT", 5) == 4
    @test randstrobe_score(dna"ATGC", dna"GTGT", 7) == 6

    @testset "get_strobe_2_mer" begin
        @test get_strobe_2_mer(dna"ATCTCTGTTT") == dna"AT--CT----"
        @test get_strobe_2_mer(test_seq) == dna"ATGC----"
        
        @test get_strobe_2_mer(dna"ATCTCTGTTT"; withGap = false) == dna"ATCT"
        @test get_strobe_2_mer(test_seq; withGap = false) == dna"ATGC"
    end
    
    @testset "ungapped_strobe_2_mer_count" begin
        counts = ungapped_strobe_2_mer_count(test_seq,s=1,w_min=2,w_max=4)
        @test round(mean(counts), digits = 4) == 0.0195
        @test counts[4] == 2
        @test counts[5] == counts[12] == counts[15] == 1
    end
end