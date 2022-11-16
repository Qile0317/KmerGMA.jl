##
#reader = open(FASTA.Reader, "C:/Users/lu_41/Desktop/Sofo Prok/VicPac32.fna")
#path = "C:/Users/lu_41/Desktop/Sofo Prok/VicPac32.fna"
#firstr = first(reader)
#f = FASTA.sequence(firstr)
#q = dna"agtatcactaattatcagagaaatgcaaatcaaaactac"
#query = ExactSearchQuery(q)
#bseq = dna"CCCCCC"^10

## I can also make it return a dictionary.
#ig also its possible to make it write into a file
"""
    firstMatch(reader,query)

Scans through a FASTA,Reader object to find the FIRST occurence of the query(dna longsequence) and prints the results to the REPL.

This function is only really needed to quickly see if there is alot of matches.
"""
function firstMatch(reader::FASTX.FASTA.Reader{}, query::LongSequence{DNAAlphabet{4}})
    query = ExactSearchQuery(query)
    for record in reader
        find = findfirst(query, FASTA.sequence(LongSequence{DNAAlphabet{4}}, record))
        if !isnothing(find)
            println(find," ", FASTA.identifier(record))
        end
    end
end

export FirstMatch

function FindAll(q::ExactSearchQuery{typeof(isequal), LongSequence{DNAAlphabet{4}}}, seq::LongSequence{DNAAlphabet{4}}, answer::Vector{UnitRange} = UnitRange[])
    start = 1
    rg = findfirst(q, view(seq, start: length(seq)))
    while !isnothing(rg)
        push!(answer, start-1+first(rg): start-1+last(rg))
        start += last(rg)
        rg = findfirst(q, view(seq, start: length(seq)))
    end
    if answer == UnitRange[]
        return nothing
    else
        return answer
    end
end

#overlap
function FindAllOverlap(q::ExactSearchQuery{typeof(isequal), LongSequence{DNAAlphabet{4}}}, seq::LongSequence{DNAAlphabet{4}}, answer::Vector{UnitRange})
    start = 1
    rg = findfirst(q, view(seq, start: length(seq)))
    while !isnothing(rg)
        push!(answer, start-1+first(rg): start-1+last(rg))
        start += first(rg)
        rg = findfirst(q, view(seq, start: length(seq)))
    end
    if answer == UnitRange[]
        return nothing
    else
        return answer
    end
end

"""
    exactMatch(query::LongSequence{DNAAlphabet{4}},
               seq::LongSequence{DNAAlphabet{4}},
               overlap::Bool = true)

Finds all exact matches to a query sequence(dna longsequence) in the given genome assembly as a reader object(seq) or single sequence

overlap is a boolean argument and is true by default

Returns a dictionary of the identifiers of individual records it found matches in and the match locations.

The algorithm is simply based on the Biosequences findfirst() function and runs quite fast through entire genomes.
"""
function exactMatch(query::LongSequence{DNAAlphabet{4}}, seq::LongSequence{DNAAlphabet{4}}; overlap::Bool = true)
    q = ExactSearchQuery(query)
    answer = UnitRange[]
    if overlap==true
        FindAllOverlap(q, seq, answer)
    elseif overlap == false
        FindAll(q,seq, answer)
    end
end

function exactMatch(query::LongSequence{DNAAlphabet{4}}, Reader::FASTX.FASTA.Reader{}; overlap::Bool = true)
    q = ExactSearchQuery(query)
    identify = Dict{String,Vector{UnitRange{Int64}}}()
    for record in Reader
        seq = FASTA.sequence(record)
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

##

#FASTX.FASTA.seqlen(firstr) # 0.000001 seconds
#is much faster than other methods

"""
    cflength(reader)

function to record and store cumulative lengths of the BEGINNING of each record in a dictionary
for example: for the reader

    >firstseq
    ATGC
    >secondseq
    AT

the dictionary returned would be:

    "firstseq"  => 4
    "secondseq" => 6

Not nessecarily in that order.
"""
function cflength(reader::FASTX.FASTA.Reader{})
    lengthmap = Dict{String,Int64}()
    clength = 0
    prev = 0
    for record in reader
        identifier = FASTA.identifier(record)
        clength += prev
        prev = FASTX.FASTA.seqlen(record)
        lengthmap[identifier] = clength
    end
    return lengthmap
end

export cflength

function PltQueryMatches(matches::Dict{String, Vector{UnitRange{Int64}}}, rlength::Int64, clengths::Dict{String,Int64})
    x = Int64[]
    for match in matches
        currlen = clengths[first(match)]
        c = matches[first(match)]
        for ur in c
            f = first(ur) + currlen
            push!(x,f)
        end
    end
    push!(x,rlength) #for a slightly more elegant graph i could maybe try find a way to extend the x axis to rlength
    y = fill(1,rlength)
    scatter(x,y, title = "query match locations", label = "first position of match")
    xlabel!("matches along genome")
end

"""
    PlotQueryMatches(matches, reader)

matches is the dictionary returned by exactMatch() from the original reader

plot the locations along the genome where an exact match was found.
"""
function PlotQueryMatches(matches::Dict{String, Vector{UnitRange{Int64}}}, reader::FASTX.FASTA.Reader)
    cl = cflength(reader)
    rlen = readerNTs(reader)
    PltQueryMatches(matches,rlen,cl)
end

export PlotQueryMatches

"""
    MatchPlot(q::LongSequence{DNAAlphabet{4}},
    Reader::FASTX.FASTA.Reader{};
    overlap::Bool)
    
function that does everything by finding matches AND plotting. doesnt work atm.
"""
function MatchPlot(q::LongSequence{DNAAlphabet{4}}, LongSequence{DNAAlphabet{4}}}, Reader::FASTX.FASTA.Reader{}; overlap::Bool)
    matches = exactMatch(q,reader,overlap)
    SlowerPQM(matches,reader)
    return matches
end

export MatchPlots

#MatchPlot(query, reader, true)
