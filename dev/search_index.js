var documenterSearchIndex = {"docs":
[{"location":"exactMatch/#Finding-exact-matches-of-a-query-in-a-fasta_file","page":"Exact gene finding","title":"Finding exact matches of a query in a fasta_file","text":"","category":"section"},{"location":"exactMatch/","page":"Exact gene finding","title":"Exact gene finding","text":"KmerGMA.firstMatch\nKmerGMA.exactMatch","category":"page"},{"location":"exactMatch/#KmerGMA.firstMatch","page":"Exact gene finding","title":"KmerGMA.firstMatch","text":"firstMatch(readerFASTX.FASTA.Reader, query::LongSequence{DNAAlphabet{4}})\n\nScans through a FASTA,Reader object to find the FIRST occurence of the query(dna longsequence) and prints the results to the REPL.\n\nThis function is only really needed to quickly see if there is alot of matches.\n\n\n\n\n\n","category":"function"},{"location":"exactMatch/#KmerGMA.exactMatch","page":"Exact gene finding","title":"KmerGMA.exactMatch","text":"exactMatch(query, subject_seq, overlap::Bool = true)\n\nFinds all exact matches to a query sequence(dna longsequence) in the given genome assembly as a reader object(seq) or single sequence\n\nquery can be a FASTA record, a substring of a dna sequence or a dna longsequence.\n\nseq can be a fasta record, dna (sub)sequence, or fasta READER.\n\noverlap is a boolean argument and is true by default\n\nReturns a dictionary of the identifiers of individual records it found matches in and the match locations.\n\nThe algorithm is simply based on the Biosequences findfirst() function and runs quite fast through entire genomes.\n\n\n\n\n\n","category":"function"},{"location":"findGenes/#Homology-searching","page":"Homology searching","title":"Homology searching","text":"","category":"section"},{"location":"findGenes/","page":"Homology searching","title":"Homology searching","text":"Here is the main function of the package. Based on a novel kmer-based algorithm to reduce runtime over a whole genome to be O(#bp) on average.","category":"page"},{"location":"findGenes/","page":"Homology searching","title":"Homology searching","text":"note: Note\nThe documentation is unfinished and may be somewhat unclear. For more details on how the algorithm functions, some pre-print of the homology searching approach will be uploaded soon.","category":"page"},{"location":"findGenes/","page":"Homology searching","title":"Homology searching","text":"The algorithm intakes a set of similar reference sequences - for example a collection of V-genes - and iterates over a genome assembly(ies) to find genes that are homologous to the references. ","category":"page"},{"location":"findGenes/","page":"Homology searching","title":"Homology searching","text":"KmerGMA.findGenes","category":"page"},{"location":"findGenes/#KmerGMA.findGenes","page":"Homology searching","title":"KmerGMA.findGenes","text":"KmerGMA.findGenes(;\n    genome_path::String,\n    ref_path::String,\n    k::Int = 6,\n    KmerDistThr::Union{Int64, Float64} = 0,\n    buffer::Int64 = 50,\n    do_align::Bool = true,\n    gap_open_score::Int = -69,\n    gap_extend_score::Int = -1,\n    do_return_dists::Bool = false,\n    do_return_hit_loci::Bool = false,\n    do_return_align::Bool = false\n    verbose::Bool = true)\n\nThe main API to conduct homology searching in a genome, using a kmer-based sequence similarity metric, against a a reference sequence set. For example, the set of all germline V genes of a mammal. Returns the approximate matches as FASTA record vector WITHIN A VECTOR of length 1 in the default configuration of the parameters. The descriptions of the record contain information about the match.  The format of the description appears as so:\n\n\"Identifier | dist = a | MatchPos = b:c | GenomePos = e | Len = f\"\n\nWhere Identifier is the contig ID of the genome where the current hit was found. dist is the kmer similarity, and MatchPos is the unitrange for the match in the contig. GenomePos is the cumulative nucleotides that have been iterated over until the current contig, and Len is simply the length of the hit length.\n\n...\n\nArguments\n\ngenome_path: a string that should be the path of the genome to conduct homology searching on.\nref_path: the file location of a fasta file containing the reference sequence(s). The references should be very similar in length.\n\nOptional Arguments\n\nk::Int64 = 6: the kmer length to use for approximate matching. Generally it should probably remain between 5 - 10 and probably does not have a profound impact on how the matches are found. Check out the pre-print for more information\nKmerDistThr::Int64 = 0: the Kmer Distance distance threshold for sequence matches to the query reference sequence set. Out of the context of the algorithm, lower values mean matches have to be more similar to the references. If left as 0, it is automatically computed. Once again the pre-print has information on this argument.\nbuffer::Int64 = 50: the amount of nucleotides left and right of a matched sequence that should be added to the returned fasta sequences, as KmerGMA is a heuristic. If do_align is set to true the buffer will be included for a more accurate alignment.\nverbose::Bool = true Indicates whether to show info in the REPL about the the progress of the processing\ndo_align: Whether to align the hits+buffer region to the consensus sequence of the references. Highly recommended to keep as true\ngap_open_score::Int = -69: The gap_open scoring paramater for the BioAlignments semi global pairwise aligner function. The number should ideally be kept low to heavily penalize gap extensions for most use-cases.\ngap_extend_score::Int = -1: The gap_extend scoring parameter for the BioAlignments semi global pairwise aligner function. Can be kept to the default -1 as long as the gap_open score is low.\ndo_return_dists: boolean to indicate whether the kmer distances along every window along the genome should be returned in a vector. (intensive memory consumption when genomes are large)\ndo_return_hit_loci: if true, will return an additional vector of the position within the genomic sequences of each hit, corresponding to the index in the hit vector.\ndo_return_align: if true, will return an additional vector of alignment object of each hit to the consensus reference sequence.\nKmerDist_threshold_buffer::Real = 8.0: a value to determine approximately how much kmer distance a hit should be lower than any random non-hit sequence. Keep in mind that kmerdistance approximates edit distance for mutations better than indels.\n\n...\n\nThe last three arguments would add term to the output. The output vector would incorporate the respective vectors in the same order of priority if any of the parameters are true.\n\nNote: Playing with the KmerDistThr argument could return more or less matches every time. \n\n\n\n\n\n","category":"function"},{"location":"findGenes/","page":"Homology searching","title":"Homology searching","text":"An alternative version of KmerGMA.findGenes which clusters the reference sequence set into similar subsets is also avaliable, though its running speed would be \"the number of subsets\" times slower. However, the following function is recommended for when speed is less important than accuracy.","category":"page"},{"location":"findGenes/","page":"Homology searching","title":"Homology searching","text":"KmerGMA.findGenes_cluster_mode","category":"page"},{"location":"findGenes/#KmerGMA.findGenes_cluster_mode","page":"Homology searching","title":"KmerGMA.findGenes_cluster_mode","text":"KmerGMA.findGenes_cluster_mode(;\n    genome_path::String,\n    ref_path::String,\n    cluster_cutoffs = [7,12,20,25],\n    k::Int = 6, \n    KmerDistThrs = Float64[0], \n    buff::Int64 = 50,\n    do_align::Bool = true,\n    gap_open_score::Int = -200,\n    gap_extend_score::Int = -1,\n    do_return_align::Bool = false,\n    do_return_dists::Bool = false,\n    do_return_hit_loci::Bool = false,\n    verbose::Bool = true,\n    kmerDist_threshold_buffer::Real = 7\n    )\n\nA slower (O(mn)) alternative to KmerGMA.findGenes to conduct homology searching in a genome, against a a reference sequence set, using a kmer-based sequence similarity metric. For example, the set of all germline V genes of a mammal.\n\nReturns the approximate matches as FASTA record vector WITHIN A VECTOR of length 1 in the default configuration of the parameters. The descriptions of the record contain information about the match.  The format of the description appears as so:\n\n\"Identifier | dist = a | KFV = b | MatchPos = c:d | GenomePos = e | Len = f\"\n\nWhere Identifier is the contig ID of the genome where the current hit was found. dist is the kmer similarity, KFV is the reference kmer frequency vector that a hit was matched against (more info below), MatchPos is the unitrange for the match in the contig. GenomePos is the cumulative nucleotides that have been iterated over until the current contig, and Len is simply the length of the hit sequence\n\n...\n\nArguments\n\ngenome_path: a string that should be the path of the genome to conduct homology searching on.\nref_path: the file location of a fasta file containing the reference sequence(s). The references should be very similar in length.\n\nOptional Arguments\n\ncluster_cutoffs = [7,12,20,25]: Important determination of the cutoff points to construct subsets of similar reference sequences. The cutoff points refer to the kmer distance against the reference. Usually, the more cutoffs may be better but the default cutoffs have been determined to be relatively robust.\nk::Int64 = 6: the kmer length to use for approximate matching. Generally it should probably remain between 5 - 10 and probably does not have a profound impact on how the matches are found. Check out the pre-print for more information\nKmerDistThrs = [0]: the Kmer Distance distance thresholds for each sequence matche to the query reference sequence set. Out of the context of the algorithm, lower values mean matches have to be more similar to the references. If left as 0, it is automatically computed. Once again a pre-print has information on this argument.\nbuffer::Int64 = 100: the amount of nucleotides left and right of a matched sequence that should be added to the returned fasta sequences, as KmerGMA is a heuristic. If do_align is set to true the buffer will be included for a more accurate alignment.\ndo_align: Whether to align the hits+bufer region to the consensus sequence of the references. Highly recommended to keep as true\ngap_open_score::Int = -200: The gap_open scoring paramater for the BioAlignments semi global pairwise aligner function. The number should ideally be kept low to heavily penalize gap extensions for most use-cases.\ngap_extend_score::Int = -1: The gap_extend scoring parameter for the BioAlignments semi global pairwise aligner function. Can be kept to the default -1 as long as the gap_open score is low.\ndo_return_align: if true, will return an additional vector of alignment object of each hit to the consensus reference sequence.\ndo_return_dists: boolean to indicate whether the kmer distances along every window along the genome should be returned in a vector. (intensive memory consumption when genomes are large)\ndo_return_hit_loci: if true, will return an additional vector of the position within the genomic sequences of each hit, corresponding to the index in the hit vector.\nverbose::Bool = true Indicates whether to show info in the REPL about the the progress of the processing\nkmerDist_threshold_buffer::Real = 7: a value to determine approximately how much kmer distance a hit should be lower than any random non-hit sequence. Keep in mind that kmerdistance approximates edit distance for mutations better than indels. Additionally, the false positive rate would likely increase as the value goes too low. Unfortunately this value is currently only optimized for k = 6. something around 20 seem to work well for k = 5. \n\n...\n\nThe last three arguments would add terms to the output. When verbose is true the exact outputs are even stated.  The output vector would incorporate the respective vectors in the same order of priority if any of the parameters are true.\n\nThere are some more details to be added in future releases.\n\n\n\n\n\n","category":"function"},{"location":"findGenes/","page":"Homology searching","title":"Homology searching","text":"To write resulting vector of hits into a fasta file, use the following function","category":"page"},{"location":"findGenes/","page":"Homology searching","title":"Homology searching","text":"KmerGMA.write_results","category":"page"},{"location":"findGenes/#KmerGMA.write_results","page":"Homology searching","title":"KmerGMA.write_results","text":"write_results(KmerGMA_result_vec::Vector{FASTX.FASTA.Record}, file_path::String, width::Int64 = 95)\n\nWrites all FASTA Records (from the FASTX package) from a vector into a fasta file location file_path. width is an optional argument that indicates the maximum width of sequences written to the file per line.\n\n\n\n\n\n","category":"function"},{"location":"findGenes/","page":"Homology searching","title":"Homology searching","text":"There is also possible pre-processing of the set of reference sequences that can be used in kmerGMA.findGenes_cluster_mode to further improve accuracy (though usually at the cost of a bit more more speed) More on this will be expanded on.","category":"page"},{"location":"findGenes/#Experimental-Strobemer-based-homology-searching","page":"Homology searching","title":"Experimental Strobemer-based homology searching","text":"","category":"section"},{"location":"findGenes/","page":"Homology searching","title":"Homology searching","text":"note: Note\nThe following function also does homology searching with strobemers but is not user-friendly and lacking more documentation and testing. ","category":"page"},{"location":"findGenes/","page":"Homology searching","title":"Homology searching","text":"KmerGMA.Strobemer_findGenes","category":"page"},{"location":"findGenes/#KmerGMA.Strobemer_findGenes","page":"Homology searching","title":"KmerGMA.Strobemer_findGenes","text":"Strobemer_findGenes(;\n    genome_path::String,\n    ref_path::String,\n    s::Int = 2,\n    w_min::Int = 3,\n    w_max::Int = 5,\n    q::Int = 5,\n    KmerDistThr::Union{Int64, Float64} = 30,\n    buffer::Int64 = 50,\n    do_align::Bool = true,\n    align_score_thr::Int = 0\n    do_return_dists::Bool = false,\n    do_return_hit_loci::Bool = false,\n    do_return_align::Bool = false,\n    verbose::Bool = true)\n\nAn experimental homology searcher that uses Strobemers, specifically randstrobes with 2 sub-kmers.  Very unoptimized and there is not distance threshold estimation yet. More documentation is to come.\n\nThe function uses the same arguments and returns the same outputs as KmerGMA.findGenes and KmerGMA.findGenes_cluster_mode but s, w_min, w_max, q are randstrobe parameters.\n\n\n\n\n\n","category":"function"},{"location":"utilities/#Utility-functions-for-processing-sequences-and-subsequences","page":"Utilities","title":"Utility functions for processing sequences and subsequences","text":"","category":"section"},{"location":"utilities/#fasta-file-utilities","page":"Utilities","title":"fasta file utilities","text":"","category":"section"},{"location":"utilities/","page":"Utilities","title":"Utilities","text":"KmerGMA.fasta_id_to_cumulative_len_dict\nKmerGMA.getSeq","category":"page"},{"location":"utilities/#KmerGMA.fasta_id_to_cumulative_len_dict","page":"Utilities","title":"KmerGMA.fasta_id_to_cumulative_len_dict","text":"fasta_id_to_cumulative_len_dict(fasta_file_path::String)\n\nfunction to record and store cumulative lengths of the BEGINNING of each record in a dictionary, from a fasta file. for example: for the following fasta file\n\n>firstseq\nATGC\n>secondseq\nAT\n\nThe dictionary returned would be:\n\n\"firstseq\"  => 4\n\"secondseq\" => 6\n\nNot nessecarily in that order.\n\n\n\n\n\n","category":"function"},{"location":"utilities/#KmerGMA.getSeq","page":"Utilities","title":"KmerGMA.getSeq","text":"getSeq(sequence::FASTA.Record)\n\nAlias to get the dna longsequence of a fasta record\n\n\n\n\n\n","category":"function"},{"location":"utilities/#kmer-utilities","page":"Utilities","title":"kmer utilities","text":"","category":"section"},{"location":"utilities/","page":"Utilities","title":"Utilities","text":"KmerGMA.kmer_count\nKmerGMA.as_kmer\nKmerGMA.as_UInt\nKmerGMA.kmer_dist","category":"page"},{"location":"utilities/#KmerGMA.kmer_count","page":"Utilities","title":"KmerGMA.kmer_count","text":"kmer_count(str::DnaSeq, k::Int, Nt_bits::DnaBits = NUCLEOTIDE_BITS)\n\nSimple kmer counting function for a DNA BioSequence str, where k is the kmer length to count. The function returns a kmer frequency vector where each INDEX of the vector corresponds to a unique kmer. For each index, the unsigned bits of the index correspond to a binary representation of sequences where two bits encode one nucleotide. To see what kmer each index corresponds to, see the function as_kmer\n\n\n\n\n\n","category":"function"},{"location":"utilities/#KmerGMA.as_kmer","page":"Utilities","title":"KmerGMA.as_kmer","text":"as_kmer(kmer_uint::Integer, kmer_len::Int,Nt_bits::Dict{UInt, Seq} = BitNtDict)\n\nTakes a positive integer representing a 2-bit-based dna kmer, and the actual length of the kmer, and returns the BioSequences dna sequence of the kmer.  The Nt_bits argument can be ignored. Additionally, note that both input parameters will be modified to 0 or 1.\n\n\n\n\n\n","category":"function"},{"location":"utilities/#KmerGMA.as_UInt","page":"Utilities","title":"KmerGMA.as_UInt","text":"UInt(kmer_seq::LongSequence{DNAAlphabet{4}}, Nt_bits::DnaBits = NUCLEOTIDE_BITS)\n\nConvert a dna biosequence kmer_seq into an unsigned integer. Users can ignore the second argument\n\n\n\n\n\n","category":"function"},{"location":"utilities/#KmerGMA.kmer_dist","page":"Utilities","title":"KmerGMA.kmer_dist","text":"kmer_dist(seq1, seq2, k::Int, Nt_bits::DnaBits = NUCLEOTIDE_BITS)\n\nreturns the kmer distance between two biosequences seq1 and seq2, where k is the kmer length. Users can ignore the last argument.\n\n\n\n\n\n","category":"function"},{"location":"utilities/#paired-kmer-utilities","page":"Utilities","title":"paired kmer utilities","text":"","category":"section"},{"location":"utilities/","page":"Utilities","title":"Utilities","text":"The following function count all kmer pairs within a sequence, which is a feature to be further utilized in an upcoming release. Counting paired kmers could serve as a more robust alternative to regular kmers in many applications.","category":"page"},{"location":"utilities/","page":"Utilities","title":"Utilities","text":"KmerGMA.kmer_pair_count","category":"page"},{"location":"utilities/#KmerGMA.kmer_pair_count","page":"Utilities","title":"KmerGMA.kmer_pair_count","text":"kmer_pair_count(seq::DnaSeq, k::Int = 3, Nt_bits::DnaBits = NUCLEOTIDE_BITS)\n\nCounts all paired kmers in a Biosequences DNA longsequence (or a longsubseq), where each index of the resulting paired kmer frequency vector correspond to the first and second kmers appended together.\n\nFor example, the 2mer pair AT and GC in a sequence ATGC would have a count of 1 in the resulting vector of length 4^(2*2) = 1 It would be at the index in the vector equivalent to as_UInt(dna\"ATGC\")\n\n\n\n\n\n","category":"function"},{"location":"utilities/#strobemer-utilities-(experimental)","page":"Utilities","title":"strobemer utilities (experimental)","text":"","category":"section"},{"location":"utilities/","page":"Utilities","title":"Utilities","text":"Strobemers consist of two or more linked shorter k-mers, where the combination of linked k-mers is decided by a hash function. This implementation is currently being experimented with and is unoptimized, but may be implemented in the near future as a new homology searching algorithm. ","category":"page"},{"location":"utilities/","page":"Utilities","title":"Utilities","text":"The strobemers implemented are simply randstrobes with two k-mers, and in the future will be a part of a new module StrobemerGMA","category":"page"},{"location":"utilities/","page":"Utilities","title":"Utilities","text":"KmerGMA.get_strobe_2_mer\nKmerGMA.ungapped_strobe_2_mer_count","category":"page"},{"location":"utilities/#KmerGMA.get_strobe_2_mer","page":"Utilities","title":"KmerGMA.get_strobe_2_mer","text":"get_strobe_2_mer(\n    seq::LongSequence{DNAAlphabet{4}},\n    s::Int = 2,\n    w_min::Int = 3,\n    w_max::Int = 5,\n    q::Int = 5;\n    withGap::Bool = true)\n\nGet the randstrobe that consists of two kmers of the current sequence seq.  The withGap named argument indicates whether to return the strobemer with sequence gaps or just as a regular kmer with all gaps removed\n\n...\n\nStrobemer parameters\n\ns::Int = 2: length of the kmers that consists of the strobemer.\nw_min::Int = 3: the start of the window to obtain the second kmer from.\nw_max::Int = 5: the end of the window to obtain the second kmer from.\nq::Int = 5: the prime number that the randstrobe hashing function should use. \n\n...\n\nIt is recommended to keep all optional strobemer parameters as is\n\n\n\n\n\n","category":"function"},{"location":"utilities/#KmerGMA.ungapped_strobe_2_mer_count","page":"Utilities","title":"KmerGMA.ungapped_strobe_2_mer_count","text":"ungapped_strobe_2_mer_count(\n    seq::DnaSeq;\n    s::Int = 2,\n    w_min::Int = 3,\n    w_max::Int = 5,\n    q::Int = 5)\n\nGet the ungapped randstrobe (two kmers of the current sequence seq) frequency vector.  The vector is 4^2s long, where each index corresponds to the strobemer with gaps removed. So the strobemer AC--GT-- would correspond to the index corresponding to ACGT. To convert from index to the ungapped strobemer, see KmerGMA.as_kmer\n\n...\n\nStrobemer parameters\n\ns::Int = 2: length of the kmers that consists of the strobemer.\nw_min::Int = 3: the start of the window to obtain the second kmer from.\nw_max::Int = 5: the end of the window to obtain the second kmer from.\nq::Int = 5: the prime number that the randstrobe hashing function should use. \n\n...\n\nIt is currently the responsibility of the user to ensure that all parameters do not conflict with eachother \n\n\n\n\n\n","category":"function"},{"location":"utilities/#sequence-divergence-(work-in-progress)","page":"Utilities","title":"sequence divergence (work in progress)","text":"","category":"section"},{"location":"utilities/","page":"Utilities","title":"Utilities","text":"KmerGMA.mutate_seq\nKmerGMA.mutate_seq!","category":"page"},{"location":"utilities/#KmerGMA.mutate_seq","page":"Utilities","title":"KmerGMA.mutate_seq","text":"mutate_seq(seq::LongSequence{DNAAlphabet{4}}, mut_rate::Real)\n\nreturns a new biosequence where mute_rate portion of the input biosequence seq are substituted to a different base pair. The mut_rate should be between 0 and 1.\n\n\n\n\n\n","category":"function"},{"location":"utilities/#KmerGMA.mutate_seq!","page":"Utilities","title":"KmerGMA.mutate_seq!","text":"mutate_seq!(seq::LongSequence{DNAAlphabet{4}}, mut_rate::Real)\n\nrandomly subsitutes mute_rate portion of the input biosequence seq in place. \n\n\n\n\n\n","category":"function"},{"location":"#KmerGMA","page":"Home","title":"KmerGMA","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"(Image: Build Status) (Image: Codecov test coverage) (Image: MIT license) (Image: Latest Release) (Image: Documentation)","category":"page"},{"location":"","page":"Home","title":"Home","text":"A k-mer based approach for homolgy searching","category":"page"},{"location":"","page":"Home","title":"Home","text":"A package for finding homologues of a reference gene family/sequence using kmer manipulation with a single scan. Also includes utilities such as an algorithm to find all exact matches of a reference sequence to a genome. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"The approach was developed to predict novel VDJ alelles in genomes based on sets of reference sequences, but should be just as viable for other gene homologues that arise from gene families.","category":"page"},{"location":"#Package-Features","page":"Home","title":"Package Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Homology searching of 1 query sequence or a collection of queries to a local database(s)\nFind exact gene matches in a genome from a query sequence\nMany subsequence counting utilities for kmers, paired kmers, and strobemers (https://doi.org/10.1101/gr.275648.121)","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The package is not yet registered in the pkg registry. Until then, please use:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg\nPkg.add(PackageSpec(name=\"KmerGMA\", url = \"https://github.com/Qile0317/KmerGMA.jl.git\"))","category":"page"},{"location":"","page":"Home","title":"Home","text":"If you are interested in the cutting edge of the development, please check out the master branch to try new features before release.","category":"page"},{"location":"#Usage-Example","page":"Home","title":"Usage Example","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To conduct homology searching of a set of sequences in a local fasta file with another query sequence fasta file, simply do:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using KmerGMA\nKmerGMA.findGenes(genome_path = \"my_sequences.fasta\", ref_path = \"my_query_sequence_family.fasta\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"Where genome_path is the file location of a fasta file containing sequences to search over, and ref_path is the file location of a fasta file containing the query sequence or a set of alike query sequences (For example V-genes). ","category":"page"},{"location":"","page":"Home","title":"Home","text":"The function defaults to returning a vector of fasta records. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"Alternatively, if accuracy is favored over speed, then its suggested to be running the following: ","category":"page"},{"location":"","page":"Home","title":"Home","text":"KmerGMA.findGenes_cluster_mode(genome_path = \"my_sequences.fasta\", ref_path = \"my_query_sequence_family.fasta\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"Where speed is sacrificed for accuracy.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The functions have many optional paramaters to optimize/adjust its performance. See the documentation for more details.","category":"page"},{"location":"#Performance","page":"Home","title":"Performance","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The current version of the homology searching findGenes function can iterate on average 40 megabases per second. So it would take about 80 seconds for the human genome. The performance of findGenes_cluster_mode slows porportionally to the number of reference sequence clusters, so for 5 clusters it would be 40/5 = 8 megabases/second. ","category":"page"},{"location":"#Contributions","page":"Home","title":"Contributions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The work was initially begun at Karolinska Institutet, as a side project to the in-progress project Discovery of Novel Germline Immunoglobulin alleles in which 2 approaches were utilized to expand camelid V(D)J databases. Thanks to @murrellb for massive support. More information is found at https://github.com/Qile0317/SoFoCompBio22","category":"page"}]
}
