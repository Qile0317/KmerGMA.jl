var documenterSearchIndex = {"docs":
[{"location":"exactMatch/#Finding-exact-match-of-a-query-in-a-genome-(fasta-file)","page":"Exact gene finding","title":"Finding exact match of a query in a genome (fasta file)","text":"","category":"section"},{"location":"exactMatch/","page":"Exact gene finding","title":"Exact gene finding","text":"The main function is simply based off of BioSequences' findfirst function.","category":"page"},{"location":"exactMatch/","page":"Exact gene finding","title":"Exact gene finding","text":"KmerGMA.exactMatch\nKmerGMA.firstMatch","category":"page"},{"location":"exactMatch/#KmerGMA.exactMatch","page":"Exact gene finding","title":"KmerGMA.exactMatch","text":"exactMatch(query,\n           seq,\n           overlap::Bool = true)\n\nFinds all exact matches to a query sequence(dna longsequence) in the given genome assembly as a reader object(seq) or single sequence\n\nquery can be a FASTA record, a substring of a dna sequence or a dna longsequence.\n\nseq can be a fasta record, dna (sub)sequence, or fasta READER.\n\noverlap is a boolean argument and is true by default\n\nReturns a dictionary of the identifiers of individual records it found matches in and the match locations.\n\nThe algorithm is simply based on the Biosequences findfirst() function and runs quite fast through entire genomes.\n\n\n\n\n\n","category":"function"},{"location":"exactMatch/#KmerGMA.firstMatch","page":"Exact gene finding","title":"KmerGMA.firstMatch","text":"firstMatch(readerFASTX.FASTA.Reader, query::LongSequence{DNAAlphabet{4}})\n\nScans through a FASTA,Reader object to find the FIRST occurence of the query(dna longsequence) and prints the results to the REPL.\n\nThis function is only really needed to quickly see if there is alot of matches.\n\n\n\n\n\n","category":"function"},{"location":"findGenes/#Homology-searching","page":"Homology searching","title":"Homology searching","text":"","category":"section"},{"location":"findGenes/","page":"Homology searching","title":"Homology searching","text":"The main function of the package. Based on a novel kmer-based algorithm to reduce runtime over a whole genome to be O(#bp).","category":"page"},{"location":"findGenes/","page":"Homology searching","title":"Homology searching","text":"The algorithm intakes a set of similar reference sequences - for example a collection of V-genes - and iterates over a genome assembly(ies) to find genes that are homologous to the references. ","category":"page"},{"location":"findGenes/","page":"Homology searching","title":"Homology searching","text":"KmerGMA.findGenes\nKmerGMA.write_results","category":"page"},{"location":"findGenes/#KmerGMA.findGenes","page":"Homology searching","title":"KmerGMA.findGenes","text":"KmerGMA.findGenes(;\n    genome_path::String,\n    ref_path::String,\n    do_cluster::Bool = false, # to be implemented\n    k::Int = 6,\n    KmerDistThr::Union{Int64, Float64} = 0,\n    buffer::Int64 = 50,\n    do_align::Bool = true,\n    do_return_dists::Bool = false,\n    do_return_hit_loci::Bool = false,\n    do_return_align::Bool = false)\n\nThe main API to conduct homolohy searching in a genome, using a kmer-based sequence similarity metric, against a a reference sequence set. For example, the set of all germline V genes of a mammal. Returns the approximate matches as FASTA record vector WITHIN A VECTOR of length 1 in the default configuration of the parameters. The descriptions of the record contain information about the match.  The format of the description appears as so:\n\n\"Identifier | dist = a | MatchPos = b:c | GenomePos = e\"\n\nWhere Identifier is the contig ID of the genome where the current hit was found. dist is the kmer similarity, and MatchPos is the unitrange for the match in the contig. GenomePos is the cumulative nucleotides that have been iterated over until the current contig.\n\n...\n\nArguments\n\ngenome_path: a string that should be the path of the genome to conduct homology searching on.\nref_path: the file location of a fasta file containing the reference sequence(s). The references should be very similar in length.\n\nOptional Arguments\n\ndo_cluster::Bool = false: work in progress, has no functionality at the moment\nk::Int64 = 6: the kmer length to use for approximate matching. Generally it should probably remain between 5 - 10 and probably does not have a profound impact on how the matches are found. Check out the pre-print for more information\nKmerDistThr::Int64 = 0: the Kmer Distance distance threshold for sequence matches to the query reference sequence set. Out of the context of the algorithm, lower values mean matches have to be more similar to the references. If left as 0, it is automatically computed. Once again the pre-print has information on this argument.\nbuffer::Int64 = 50: the amount of nucleotides left and right of a matched sequence that should be added to the returned fasta sequences, as KmerGMA is a heuristic\ndo_align: Whether to align the hits+bufer region to the consensus sequence of the references. Highly recommended to keep as true\ndo_return_dists: dangerous boolean to indicate whether the kmer distances along every window alon ghte genome should be returned in a vector\ndo_return_hit_loci: if true, will return an additional vector of the position within the genomic sequences of each hit, corresponding to the index in the hit vector.\ndo_return_align: if true, will return an additional vector of alignment object of each hit to the consensus reference sequence.\n\n...\n\nThe last three arguments would add term to the output. The output vector would incorporate the respective vectors in the same order of priority if any of the parameters are true.\n\nNote: Playing with the KmerDistThr argument could return more or less matches every time. \n\n\n\n\n\n","category":"function"},{"location":"findGenes/#KmerGMA.write_results","page":"Homology searching","title":"KmerGMA.write_results","text":"write_results(KmerGMA_result_vec::Vector{FASTX.FASTA.Record}, file_path::String, width::Int64 = 95)\n\nWrites all FASTA Records (from the FASTX package) from a vector into a fasta file location file_path. width is an optional argument that indicates the maximum width of sequences written to the file per line.\n\n\n\n\n\n","category":"function"},{"location":"findGenes/","page":"Homology searching","title":"Homology searching","text":"note: Note\nThe documentation is unfinished and may be somewhat unclear. For more details on how the algorithm functions, some pre-print of the homology searching approach will be uploaded soon.","category":"page"},{"location":"#img-src\"KmerGMA.jl.png\"-width\"30%\"-align\"right\"-/-KmerGMA","page":"Home","title":"<img src=\"KmerGMA.jl.png\" width=\"30%\" align=\"right\" /> KmerGMA","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"(Image: Build Status) (Image: Codecov test coverage) (Image: MIT license) (Image: Latest Release) (Image: Documentation)","category":"page"},{"location":"","page":"Home","title":"Home","text":"A k-mer based approach for homolgy searching","category":"page"},{"location":"","page":"Home","title":"Home","text":"A package for finding homologues of a reference gene family/sequence using kmer manipulation with a single scan. Also includes utilities such as an algorithm to find all exact matches of a reference sequence to a genome. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"The approach was developed to predict novel VDJ alelles in genomes based on sets of reference sequences, but should be just as viable for other gene homologues that arise from gene families.","category":"page"},{"location":"#Package-Features","page":"Home","title":"Package Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Homology searching of 1 query sequence or a collection of queries to a local database(s)\nFind exact gene matches in a genome from a query sequence","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The package is not yet registered in the pkg registry. Until then, please use:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg\nPkg.add(PackageSpec(name=\"KmerGMA\", url = \"https://github.com/Qile0317/KmerGMA.jl.git\"))","category":"page"},{"location":"","page":"Home","title":"Home","text":"If you are interested in the cutting edge of the development, please check out the master branch to try new features before release.","category":"page"},{"location":"#Usage-Example","page":"Home","title":"Usage Example","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To conduct homology searching of a set of sequences in a local fasta file with another query sequence fasta file, simply do:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using KmerGMA\nKmerGMA.findGenes(genome_path = \"my_sequences.fasta\", ref_path = \"my_query_sequence_family.fasta\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"Where genome_path is the file location of a fasta file containing sequences to search over, and ref_path is the file location of a fasta file containing the query sequence or a set of alike query sequences (For example V-genes). ","category":"page"},{"location":"","page":"Home","title":"Home","text":"The function defaults to returning a vector of fasta records. The findGenes function has many optional paramaters to optimize/adjust its performance. See the documentation for more details.","category":"page"},{"location":"#Performance","page":"Home","title":"Performance","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The current version of the homology searching findGenes function can iterate on average 40 megabases per second. So it would take about 80 seconds for the human genome.","category":"page"},{"location":"#Contributions","page":"Home","title":"Contributions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The work was initially begun at Karolinska Institutet, as a side project to the in-progress project Discovery of Novel Germline Immunoglobulin alleles in which 2 approaches were utilized to expand camelid V(D)J databases. Thanks to @murrellb for massive support. More information is found at https://github.com/Qile0317/SoFoCompBio22","category":"page"}]
}
