push!(LOAD_PATH,"../src/")
using EDKit
using Documenter
makedocs(
         sitename = "EDKit.jl",
         modules  = [EDKit],
         pages=[
                "Home" => "index.md"
               ])
deploydocs(;
    repo="github.com/jayren3996/EDKit.jl.git",
)