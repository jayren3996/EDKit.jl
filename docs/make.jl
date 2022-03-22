using Documenter, EDKit
push!(LOAD_PATH,"../src/")

makedocs(
    sitename="EDKit.jl",
    pages = [
        "Introduction" => "index.md",
        "Examples" => "expl.md",
        "Basis" => "basis.md",
        "Operator" => "opts.md"
    ]
)

