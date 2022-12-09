var documenterSearchIndex = {"docs":
[{"location":"exactMatch/#Finding-exact-match-of-a-query-in-a-genome-(fasta-file)","page":"Exact gene finding","title":"Finding exact match of a query in a genome (fasta file)","text":"","category":"section"},{"location":"exactMatch/","page":"Exact gene finding","title":"Exact gene finding","text":"The main function is simply based off of BioSequences' findfirst function.","category":"page"},{"location":"exactMatch/","page":"Exact gene finding","title":"Exact gene finding","text":"KmerGMA.exactMatch\nKmerGMA.firstMatch","category":"page"},{"location":"exactMatch/#KmerGMA.exactMatch","page":"Exact gene finding","title":"KmerGMA.exactMatch","text":"exactMatch(query,\n           seq,\n           overlap::Bool = true)\n\nFinds all exact matches to a query sequence(dna longsequence) in the given genome assembly as a reader object(seq) or single sequence\n\nquery can be a FASTA record, a substring of a dna sequence or a dna longsequence.\n\nseq can be a fasta record, dna (sub)sequence, or fasta READER.\n\noverlap is a boolean argument and is true by default\n\nReturns a dictionary of the identifiers of individual records it found matches in and the match locations.\n\nThe algorithm is simply based on the Biosequences findfirst() function and runs quite fast through entire genomes.\n\n\n\n\n\n","category":"function"},{"location":"exactMatch/#KmerGMA.firstMatch","page":"Exact gene finding","title":"KmerGMA.firstMatch","text":"firstMatch(readerFASTX.FASTA.Reader, query::LongSequence{DNAAlphabet{4}})\n\nScans through a FASTA,Reader object to find the FIRST occurence of the query(dna longsequence) and prints the results to the REPL.\n\nThis function is only really needed to quickly see if there is alot of matches.\n\n\n\n\n\n","category":"function"},{"location":"findGenes/#Homology-searching","page":"Homology searching","title":"Homology searching","text":"","category":"section"},{"location":"findGenes/","page":"Homology searching","title":"Homology searching","text":"The main function of the package. Based on a novel kmer-based algorithm to reduce runtime over a whole genome to be O(#bp).","category":"page"},{"location":"findGenes/","page":"Homology searching","title":"Homology searching","text":"The algorithm intakes a set of similar reference sequences - for example a collection of V-genes - and iterates over a genome assembly(ies) to find genes that are homologous to the references. ","category":"page"},{"location":"findGenes/","page":"Homology searching","title":"Homology searching","text":"KmerGMA.findGenes","category":"page"},{"location":"findGenes/#KmerGMA.findGenes","page":"Homology searching","title":"KmerGMA.findGenes","text":"  findGenes(; genome,\n              ref,\n              hasN::Bool = true,\n              mode::String = \"return\",\n              fileloc::String = \"noPath\",\n              k::Int64 = 6,\n              windowsize::Int64 = 0,\n              thr::Int64 = 0,\n              buffer::Int64 = 50,\n              resultVec::Vector{FASTA.Record} = FASTA.Record[])\n\nThe main API to find all gene matches in a genome from a reference sequence/reference sequence set. Returns the approximate matches as FASTA records. The descriptions of the record contain information about the match.  The format of the description appears as so:\n\n\"Identifier | SED = a | Pos = b:c | thr = d | buffer = e\"\n\nWhere SED indicates how similar the match is to the references, and the lower, the more similar in terms of kmer distance.\n\n...\n\nArguments\n\ngenome: Either a String or a Vector of strings, where each string is the file location of a .fasta file containing the genome.\nref: the file location of a fasta file containing the reference sequence(s). The references should be very similar in length. OR a pre-generated reference kmer frequency dictionary.\n\nOptional Arguments\n\nhasN::Bool = true: Whether both the fasta files from ref and the query genome contains any undefined nucleotide N. If false it runs an alternate algorithm that is alot faster than the one accounting for N nts. \nmode = \"return\": what the function should return the gene matches in. There are three modes, return, print, write. \n\nReturn mode returns a vector of FASTA records. Print mode simply prints the fasta sequences in the REPL.  Write mode writes the fasta sequences to a fasta file, of which the location has to be given in the fileloc argument\n\nfileloc::String = \"noPath\": If the mode argument is \"write\", this should be the location of a fasta file that the results should be written into.\nk::Int64 = 6: the kmer length to use for approximate matching. Generally it should probably remain between 5 - 10 and probably does not have a profound impact on how the matches are found. Check out the pre-print for more information\nwindowsize::Int64 = 0: the size of the returned matches. It should be the average length of all reference sequences and is computed automatically if the argument is left as 0.\nthr::Int64 = 0: the squared euclidian distance threshold for sequences. Out of the context of the algorithm, lower values mean matches have to be more similar to the references. If left as 0, it is automatically computed. Once again the pre-print has information on this argument.\nbuffer::Int64 = 0: the amount of nucleotides left and right of a matched sequence that should be added to the returned fasta sequences\n\n...\n\nIt is recommended to leave all other optional arguments as is, especially the buffer = 50 argument. The purpose of the algorithm is to find approximate matches that can then be BLASTed and aligned.\n\nHowever, playing with the thr argument could return more or less matches every time. \n\n\n\n\n\n","category":"function"},{"location":"findGenes/","page":"Homology searching","title":"Homology searching","text":"As mentioned in the docstring, if there surely are no \"N\" nucleotides present, setting HasN to false runs an alternate, much faster algorithm.","category":"page"},{"location":"findGenes/","page":"Homology searching","title":"Homology searching","text":"The thr argument can be toyed around with, as increasing the kmer distance threshold doesn't nessecarily always increase the amount of matches. The default distance is approximate against random sequences and may not always be the best metric.","category":"page"},{"location":"findGenes/","page":"Homology searching","title":"Homology searching","text":"note: Note\nThe unoptimized algorithm that accounts for N nucleotides runs in just under 5 minutes to iterate over a genome of almost 4 billion bps. More benchmarking is on the way. More information is upcoming in a pre-print.","category":"page"},{"location":"Utils/#Utilities","page":"Utilities","title":"Utilities","text":"","category":"section"},{"location":"Utils/","page":"Utilities","title":"Utilities","text":"These are a random collection of functions to do with counting sizes, record sizes, and dictionary manipulation functions. ","category":"page"},{"location":"Utils/","page":"Utilities","title":"Utilities","text":"note: Note\nThey are just here as this is a preliminary version of the package. ","category":"page"},{"location":"Utils/#Kmer-dictionary-generation","page":"Utilities","title":"Kmer dictionary generation","text":"","category":"section"},{"location":"Utils/","page":"Utilities","title":"Utilities","text":"KmerGMA.genKmers","category":"page"},{"location":"Utils/#KmerGMA.genKmers","page":"Utilities","title":"KmerGMA.genKmers","text":"genKmers(k::Int64,\n         Dictionary::Bool = true;\n         withN::Bool = false,\n         cumulative::Bool = false)\n\nsystematic Kmer generation into dictionary：It matches each kmer to a unique index.\n\nthe withN argument can be specified to decide whether dna\"N\" should be included as a nucleotide in the resulting kmers.\n\nthe cumulative argument shouldnt ever be needed but it generates a kmer Dictionary with all the kmer sizes up to k.\n\n\n\n\n\n","category":"function"},{"location":"Utils/#Rudimentary-kmer-counting,-very-unoptimized","page":"Utilities","title":"Rudimentary kmer counting, very unoptimized","text":"","category":"section"},{"location":"Utils/","page":"Utilities","title":"Utilities","text":"KmerGMA.kmerFreq","category":"page"},{"location":"Utils/#KmerGMA.kmerFreq","page":"Utilities","title":"KmerGMA.kmerFreq","text":"kmerFreq(k::Int64,\n     seq::LongSequence{DNAAlphabet{4}},\n     KD::Dict{LongSequence{DNAAlphabet{4}};\n     returnDict::Bool = false)\n\nComputes the kmer frequency vector of a sequence from scratch by iteration (O(n)). It is actually a slightly slower version that the main GMA does not use.\n\nKD is the kmer Dictionary and an optional argument. If left empty, will be automatically generated WITH N's.\n\nNote that it includes overlap(it can be tested whether it helps. but intuitively the difference shouldnt be too large).\n\n\n\n\n\n","category":"function"},{"location":"Utils/#Other-utilities","page":"Utilities","title":"Other utilities","text":"","category":"section"},{"location":"Utils/","page":"Utilities","title":"Utilities","text":"KmerGMA.avgRecLen\nKmerGMA.kfv\nKmerGMA.percentN\nKmerGMA.readerNTs\nKmerGMA.readerlens\nKmerGMA.flipDict","category":"page"},{"location":"Utils/#KmerGMA.avgRecLen","page":"Utilities","title":"KmerGMA.avgRecLen","text":"avgRecLen(reader;\n          rnd::Bool = true)\n\nget the rounded average length of every record in a FASTA.Reader object. (it can also be a string indicating the path)\n\nrnd is an optional argument for whether the result should be rounded to an Interger.\n\n\n\n\n\n","category":"function"},{"location":"Utils/#KmerGMA.kfv","page":"Utilities","title":"KmerGMA.kfv","text":"kfv(maindict::Dict{LongSequence{DNAAlphabet{4}}, Float64},\n    KD::Dict{LongSequence{DNAAlphabet{4}}, Int64})\n\nConverts a dictionary of kmers and COUNTS (maindict) to a kmer frequency VECTOR.\n\nKD is the kmer dictionary and must be based on the the same kmer length as the count dictionary.\n\n\n\n\n\n","category":"function"},{"location":"Utils/#KmerGMA.percentN","page":"Utilities","title":"KmerGMA.percentN","text":"percentN(seq)\n\nSimple function to count the percentage of N nucleotides in a longsequence, record, or reader object.\n\n\n\n\n\n","category":"function"},{"location":"Utils/#KmerGMA.readerNTs","page":"Utilities","title":"KmerGMA.readerNTs","text":"readerNTs(reader::FASTX.FASTA.Reader)\n\nfunction to count the number of nts in a fasta reader object (length)\n\n\n\n\n\n","category":"function"},{"location":"Utils/#KmerGMA.readerlens","page":"Utilities","title":"KmerGMA.readerlens","text":"readerlens(reader::FASTX.FASTA.Reader)\n\nfunction to just look at the individual record lengths of a reader object and returns a vector of lengths\n\n\n\n\n\n","category":"function"},{"location":"Utils/#KmerGMA.flipDict","page":"Utilities","title":"KmerGMA.flipDict","text":"flipDict(dict::Dict{LongSequence{DNAAlphabet{4}},Int64})\n\nfunction to return a NEW dictionary that has the key value pairs flipped compared to the input.\n\n\n\n\n\n","category":"function"},{"location":"#KmerGMA.jl","page":"Home","title":"KmerGMA.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"(Image: Build Status) (Image: MIT license) (Image: Latest Release) (Image: Documentation)","category":"page"},{"location":"","page":"Home","title":"Home","text":"A k-mer based approach for homolgy searching","category":"page"},{"location":"","page":"Home","title":"Home","text":"A work-in-progress package for finding homologues of a reference gene family/sequence using kmer manipulation swiftly. Also includes utilities such as an algorithmic to find all exact matches of a reference sequence to a genome. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"note: Note\nThe package is in a very preliminary stage and can be heavily optimized. For potential developers, a pre-print of the main kmer-based algorithm will be available very soon. ","category":"page"},{"location":"#Package-Features","page":"Home","title":"Package Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Homology searching of 1 query sequence or a collection of queries to a local database(s)\nFind exact gene matches in a genome from a query sequence\nUtilities for kmer-tricks","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The package is not yet registered in the pkg registry. Until then, please use:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg\nPkg.add(PackageSpec(name=\"KmerGMA\", url = \"https://github.com/Qile0317/KmerGMA.jl.git\"))","category":"page"},{"location":"","page":"Home","title":"Home","text":"If you are interested in the cutting edge of the development, please check out the master branch to try new features before release.","category":"page"},{"location":"#Performance","page":"Home","title":"Performance","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The current version of the approximate gene matching algorithm in the package can iterate over a collection of sequences of around 4 billion basepairs in 5 minutes. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"There also is a version of the algorithm that runs approximately 3 times faster if the user knows that no sequence imput in the algorithm contains the undefined nucleotide N.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Current benchmarking indicates that optimziations can improve this by a full minute in subsequent updates. Proposed optimizations has the theoretical potential to run in under a minute. More details will be shown very soon in a pre-print. ","category":"page"}]
}
