push!(LOAD_PATH,"../src/")
using KmerGMA
using Documenter
makedocs(
         sitename = "KmerGMA.jl",
         modules  = [KmerGMA],
         pages=[
                "Home" => "index.md"
               ])
deploydocs(;
    repo="github.com/Qile0317/KmerGMA.jl",
)
