using Documenter
using EDKit
using ITensors
using ITensorMPS
using LinearAlgebra
using SparseArrays

DocMeta.setdocmeta!(
    EDKit,
    :DocTestSetup,
    :(using EDKit, ITensors, ITensorMPS, LinearAlgebra, SparseArrays),
    recursive = true,
)

makedocs(
    sitename = "EDKit.jl",
    authors = "JayRen and contributors",
    modules = [EDKit],
    checkdocs = :none,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://jayren3996.github.io/EDKit.jl",
    ),
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting-started.md",
        "Manual" => [
            "Architecture" => "manual/architecture.md",
            "Bases and Sectors" => "manual/bases.md",
            "General Abelian Symmetries" => "abelian_basis.md",
            "Operators" => "manual/operators.md",
            "Maps and Symmetrizers" => "manual/maps.md",
            "Entanglement" => "manual/entanglement.md",
            "ITensor Workflows" => "manual/itensors.md",
            "Lindblad Workflows" => "manual/lindblad.md",
            "Utilities" => "manual/utilities.md",
        ],
        "Worked Examples" => [
            "Basic Workflows" => "examples/basic-workflows.md",
            "Symmetry Workflows" => "examples/symmetry-workflows.md",
            "Tensor-Network Workflows" => "examples/tensor-network-workflows.md",
            "Open-System Workflows" => "examples/open-system-workflows.md",
        ],
        "API Reference" => [
            "Overview" => "reference/api-overview.md",
            "Bases" => "reference/bases.md",
            "Operators" => "reference/operators.md",
            "Maps" => "reference/maps.md",
            "Entanglement" => "reference/entanglement.md",
            "ITensors" => "reference/itensors.md",
            "Lindblad" => "reference/lindblad.md",
            "Utilities" => "reference/utilities.md",
        ],
        "Developer Notes" => "developer.md",
    ],
)
