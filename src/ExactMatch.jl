"""
    firstMatch(readerFASTX.FASTA.Reader, query::LongSequence{DNAAlphabet{4}})

Scans through a FASTA,Reader object to find the FIRST occurence of the query(dna longsequence) and prints the results to the REPL.

This function is only really needed to quickly see if there is alot of matches.
"""
function firstMatch(reader::FASTX.FASTA.Reader, query::LongSequence{DNAAlphabet{4}})
    query = ExactSearchQuery(query)
    for record in reader
        find = findfirst(query, FASTA.sequence(LongSequence{DNAAlphabet{4}}, record))
        if !isnothing(find)
            println(find," ", FASTA.identifier(record))
        end
    end
end

export FirstMatch

function FindAll(q::ExactSearchQuery{typeof(isequal), LongSequence{DNAAlphabet{4}}},
    seq::LongSequence{DNAAlphabet{4}}, answer::Vector{UnitRange} = UnitRange[])
    start = 1
    rg = findfirst(q, view(seq, start: length(seq)))
    while !isnothing(rg)
        push!(answer, start-1+first(rg): start-1+last(rg))
        start += last(rg)
        rg = findfirst(q, view(seq, start: length(seq)))
    end
    answer == UnitRange[] ? (return nothing) : (return answer)
end

#overlap
function FindAllOverlap(q::ExactSearchQuery{typeof(isequal), LongSequence{DNAAlphabet{4}}},
    seq::LongSequence{DNAAlphabet{4}}, answer::Vector{UnitRange})
    start = 1
    rg = findfirst(q, view(seq, start: length(seq)))
    while !isnothing(rg)
        push!(answer, start-1+first(rg): start-1+last(rg))
        start += first(rg)
        rg = findfirst(q, view(seq, start: length(seq)))
    end
    answer == UnitRange[] ? (return nothing) : (return answer)
end

"""
    exactMatch(query, seq, overlap::Bool = true)

Finds all exact matches to a query sequence(dna longsequence) in the given genome assembly as a reader object(seq) or single sequence

query can be a FASTA record, a substring of a dna sequence or a dna longsequence.

seq can be a fasta record, dna (sub)sequence, or fasta READER.

overlap is a boolean argument and is true by default

Returns a dictionary of the identifiers of individual records it found matches in and the match locations.

The algorithm is simply based on the Biosequences findfirst() function and runs quite fast through entire genomes.
"""
function exactMatch(query, seq; overlap::Bool = true)
    if typeof(query) == FASTX.FASTA.Record
        query = getSeq(query)
    end
    q = ExactSearchQuery(LongSequence{DNAAlphabet{4}}(query))
    if typeof(seq) == FASTX.FASTA.Record
        seq = getSeq(seq)
    else
        seq = LongSequence{DNAAlphabet{4}}(seq)
    end
    answer = UnitRange[]
    if overlap
        FindAllOverlap(q, seq, answer)
    else
        FindAll(q,seq, answer)
    end
end

#had to get rid of some type specifications
function exactMatch(query, Reader::FASTX.FASTA.Reader{}; overlap::Bool = true)
    if typeof(query) == FASTX.FASTA.Record
        query = getSeq(query)
    end
    q = ExactSearchQuery(LongSequence{DNAAlphabet{4}}(query))
    identify = Dict{String,Vector{UnitRange{Int64}}}()
    for record in Reader
        seq = getSeq(record)
        RM = exactMatch(query,seq; overlap=overlap)
        if !isnothing(RM)
            identify[FASTA.identifier(record)] = RM
        end
    end
    if identify == Dict{String, Vector{UnitRange{Int64}}}()
        return "no match"
    else
        return identify
    end
    close(reader)
end

export exactMatch

#the algo is nothing fancy, it just uses BioSequences's exactmatch function. It'll be interesting comparing it to the GMA at SED = 0

"""
    fasta_id_to_cumulative_len_dict(fasta_file_path::String)

function to record and store cumulative lengths of the BEGINNING of each record in a dictionary, from a fasta file.
for example: for the following fasta file

    >firstseq
    ATGC
    >secondseq
    AT

The dictionary returned would be:

    "firstseq"  => 4
    "secondseq" => 6

Not nessecarily in that order.
"""
function fasta_id_to_cumulative_len_dict(fasta_file_path::String)
    lengthmap = Dict{String,Int64}()
    clength, prev = 0, 0
    open(FASTX.FASTA.Reader, fasta_file_path) do reader
        for record in reader
            identifier = FASTA.description(record)
            clength += prev
            prev = FASTX.FASTA.seqsize(record)
            lengthmap[identifier] = clength
        end
    end
    return lengthmap
end

export fasta_id_to_cumulative_len_dict
