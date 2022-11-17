#ExactMatch plotting
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

#export PlotQueryMatches

## simpleExplore plotting


"""
    refPlot(k::Int64,
    referece::Vector{Float64};
    nonzero::Bool = false)

A not so useful function to plot a kmer dictionary with some labels. One can easily customize the plot with their own plotting functions.
"""
function refPlot(k::Int64, referece::Vector{Float64}; nonzero::Bool = false)
    if nonzero == false
        plot(collect(1:1:length(referece)),referece,label=nothing)
        xlabel!(string("All Unique ",k,"-mers"))
        ylabel!("Average counts")
        title!(string("Average ",k,"-mer distribution in the reference"))
    elseif nonzero == true
        v = Int64[]
        for i in referece
            if i == 0.0
                push!(v,0)
            else
                push!(v,1)
            end
        end
        scatter(v,label=nothing)
        xlabel!(string("All Unique ",k,"-mers"))
        ylabel!("occurence")
        title!(string("all ",k,"-mer occurences in the average of the reference"))

    end
end

function refPlot(referece::Dict{LongSequence{DNAAlphabet{4}}, Float64}; nonzero::Bool = false, k::Int64 = 0)
    if k==0
        k = length(first(first(V3Ref)))
    end
    referece = dictToVec(referece)
    if nonzero == false
        plot(collect(1:1:length(referece)),referece,label=nothing)
        xlabel!(string("All Unique ",k,"-mers"))
        ylabel!("Average counts")
        title!(string("Average ",k,"-mer distribution in the reference"))
    elseif nonzero == true
        v = Int64[]
        for i in referece
            if i == 0.0
                push!(v,0)
            else
                push!(v,1)
            end
        end
        scatter(v,label=nothing)
        xlabel!(string("All Unique ",k,"-mers"))
        ylabel!("occurence")
        title!(string("all ",k,"-mer occurences in the average of the reference"))
    end
end

#export refPlot

##Plotting results (I also need to do a different type of plot for eucledian distances.)
function freqTable(k::Int64, FreqDict::Dict{Int64,Int64}, Sort::Bool = true) #I need to make the title where you can change it according to k. Ik theres an easy way to do it but i havent looked it up yet.
    if Sort == true
        bar(sort(collect(keys(FreqDict))), collect(values(FreqDict)), orientation=:vertical, label = nothing)
    elseif Sort == false
        bar(collect(keys(Freqs)), collect(values(FreqDict)), orientation=:vertical, label = nothing)
    end
    xlabel!("All Unique kmers")
    ylabel!("Kmer Counts")
    title!(string(k,"-mer Counts in sequence"))
end

#adding method for vector.
function freqTable(FreqVec::Vector{Int64})
    k = Int64(round(sqrt(sqrt(length(FreqVec)))))
    bar((FreqVec), orientation=:vertical, label = nothing)
    xlabel!("All Unique kmers")
    ylabel!("Kmer Counts")
    title!(string(k,"-mer Counts in sequence"))
end

#I also need to plot the euclidian distances and some otherstuff
# I should alos implement an optional "threshold option to just plot a line"
function eucTable(eucVec::Vector{Float64}, thr::Union{Float64,Int64} = 0)
    e = plot(eucVec, label = nothing,linewidth=0.3)
    if thr != 0
        plot(e, fill(minimum(eucVec)+thr,length(eucVec)),label = "Threshold")
    else
        plot(eucVec, label = nothing,linewidth=0.3)
    end
    xlabel!("Position in genome(bp)")
    ylabel!("Euclidian distance to reference")
end

"""
    readerlens(reader::FASTX.FASTA.Reader; plt::Bool = true)

function to just look at the individual record lengths of a reader object and optionally plot them with the plt argument.
"""
function readerlens(reader::FASTX.FASTA.Reader; plt::Bool = true)
    lenvec = Int64[]
    for record in reader
        push!(lenvec,FASTX.FASTA.seqsize(record))
    end
    if !plt
        return lenvec
    else
        Plots.bar(lenvec, title="lengths of individual records in the FASTA file", label=nothing)
        xlabel!("sequence index")
        ylabel!("sequence length")
    end
end

#export readerlens

##From when I was distracted:


#maybe kmer spectra that fills a vector. KAT.jl has cooler ways to do this.
#problem: for something like VicPac theres too many different values. Perhaps making a histogram with each value or making some CDF is better.
function kmerSpectra(k::Int64, KCDict::Dict{LongSequence{DNAAlphabet{4}},Int64}, rm::Bool = false)
    ans = fill(0,(2^((4*k-1)))-(2^((2*k)-1)))
    for key in KCDict
        ans[last(key)+1] += 1
    end
    if rm
        Ransvec = reverse(ans)
        stop = false
        for i in Ransvec
            if stop == false
                if i == 0
                    pop!(ans)
                else
                    stop = true #there should be a better way
                end
            else
                break
            end
        end
    end
    Plots.scatter(ans,label=nothing)
    xlabel!(string(k)*"-mer frequency + 1")
    ylabel!("count")
end

function oldKmerSpectra(k::Int64, KCDict::Dict{LongSequence{DNAAlphabet{4}},Int64},typ::String="bar")
    ans = Dict{Int64,Int64}()
    @showprogress "initializing frequency array..." for key in KCDict
        if haskey(ans,last(key))
            ans[last(key)] += 1
        else
            ans[last(key)] = 1
        end
    end
    ansvec = fill(0,(2^((4*k-1)))-(2^((2*k)-1))) #this is to minimize the amount of terms to put in the vector
    for key in ans
        ansvec[first(key)+1]=last(key)
    end
    Ransvec = reverse(ansvec)
    stop = false
    for i in Ransvec
        if stop == false
            if i == 0
                pop!(ansvec)
            else
                stop = true #there should be a better way
            end
        else
            break
        end
    end
    if typ == "bar"
        Plots.bar(ansvec,label=nothing)
    else
        Plots.scatter(ansvec,label=nothing)
    end
    xlabel!(string(k)*"-mer frequency + 1")
    ylabel!("count")
end
