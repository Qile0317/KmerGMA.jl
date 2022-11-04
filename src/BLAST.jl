#Optional BLAST pipeline. Based off of Ben Murrell's WebBlast.jl
"""
    flatten_to_dataframe(blast_arr)
Takes an output array of dictionaries (ie. something returned by a ```WebBLAST(query)``` call) and flattens into a DataFrame.
"""
function flatten_to_dataframe(blast_arr)
    #The fields to keep are the things in here, minus the array fields
    #println(union(vcat(collect.(keys.(blast_arr))...)))
    #println(union(vcat(collect.(keys.(blast_arr[1]["Query_hit_array"]))...)))
    #println(union(vcat(collect.(keys.(blast_arr[1]["Query_hit_array"][1]["Hit_hsps"]))...)))
    query_fields = ["Query-def", "Query_iter-num", "Query-len"]
    hit_fields = ["Hit_accession", "Hit_def", "Hit_len", "Hit_num", "Hit_id"]
    HSP_fields = ["Hsp_query-from", "Hsp_align-len", "Hsp_hseq", "Hsp_hit-from", "Hsp_score", "Hsp_positive", "Hsp_qseq", "Hsp_query-to", "Hsp_query-frame", "Hsp_midline", "Hsp_evalue", "Hsp_hit-to", "Hsp_identity", "Hsp_hit-frame", "Hsp_gaps", "Hsp_bit-score", "Hsp_num"];
    flat_fields = vcat(query_fields,hit_fields,HSP_fields);

    field_arr_of_arr = [[] for i in 1:length(flat_fields)]
    for query in blast_arr
        for hit in query["Query_hit_array"]
            for HSP in hit["Hit_hsps"]
                ind = 1
                for f in query_fields
                    push!(field_arr_of_arr[ind],get(query,f,""))
                    ind += 1
                end
                for f in hit_fields
                    push!(field_arr_of_arr[ind],get(hit,f,""))
                    ind += 1
                end
                for f in HSP_fields
                    push!(field_arr_of_arr[ind],get(HSP,f,""))
                    ind += 1
                end
            end
        end
    end

    df = DataFrame()
    for i in 1:length(flat_fields)
        df[!,flat_fields[i]] = field_arr_of_arr[i]
    end
    return df
end

"""
    eqBLAST(path::String;
            nhits::Int64 = 10, wDesc::Bool = true, df::Bool = true,
            maxWait::Union{Int64,Float64} = Inf, verb::Int64 = 1, opstring::String = "",
            XMLpath::Union{String, Nothing} = nothing, db::String = "nt", prog::String = "blastn",
            blastURL::String = "https://blast.ncbi.nlm.nih.gov/Blast.cgi?")

using WebBLAST.jl by the murrell group, it takes the filepath of a fasta file and blasts it.
Dont use yet because its extremely slow since it converts the file to 2 vectors.
I am working on making it be able to directl blast a file.

- The optional arguments do not work at the moment and it will be fixed.
"""
function eqBLAST(path::String; nhits::Int64 = 10, wDesc::Bool = true, df::Bool = true,
    maxWait::Union{Int64,Float64} = Inf, verb::Int64 = 1, opstring::String = "",
    XMLpath::Union{String, Nothing} = nothing, db::String = "nt", prog::String = "blastn",
    blastURL::String = "https://blast.ncbi.nlm.nih.gov/Blast.cgi?")
    #conversion into sequence string vector
    reader = open(FASTX.FASTA.Reader, path)
    queries = String[]
    qnames = String[]
    for record in reader
        push!(queries, string(FASTX.FASTA.sequence(record)))
        if wDesc
            push!(qnames, string(FASTX.FASTA.identifier(record)*FASTX.FASTA.description(record)))
        else
            push!(qnames, string(FASTX.FASTA.identifier(record)))
        end
    end
    #blasting
    blast_array = WebBlast.WebBLAST(queries, query_names = qnames) #, max_waits = maxWait,
    #num_hits = nhits, database = db, program = prog, verbosity = verb,
    #option_string = opstring, save_XML_path = XMLpath, Blast_URL = blastURL)
    if !df
        return blast_array
    else
        return flatten_to_dataframe(blast_array)
    end
end

export eqBLAST

#@time hits = eqBLAST(vpscan)
#names(hits)
#evals = hits[:, "Hsp_evalue"]

#@time GMA(VicPac,AlpacaV,"VicPacScan/vicpacscan.fasta")
# BoundsError: attempt to access 15625-element Vector{Float64} at index [0]
#its fine bc last record in vicpac is super short.

#using Distributed #i cant seem to get it to work atm

#phylo stuff. ALso bens package webblast allows blast within julia
#using PhyloNetworks, PhyloPlots

#interesting to note that even with different referrences that were actually valid, it had near identical performance
#this might be interesting to put into the paper.
