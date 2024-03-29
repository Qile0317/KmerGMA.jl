push!(LOAD_PATH,"../src/")
using KmerGMA
using Documenter
makedocs(
         sitename = "KmerGMA.jl",
         modules  = [KmerGMA],
         pages=[
                "Home" => "index.md"
                "Homology searching" => "findGenes.md"
                "Exact gene finding" => "exactMatch.md"
                "Utilities" => "utilities.md"
               ])
deploydocs(;
    repo="github.com/Qile0317/KmerGMA.jl.git", branch = "gh-pages", devbranch = "master"
)