using FVM_PDEsolver
using Documenter

DocMeta.setdocmeta!(FVM_PDEsolver, :DocTestSetup, :(using FVM_PDEsolver); recursive=true)

makedocs(;
    modules=[FVM_PDEsolver],
    authors="Amir Farzin",
    repo="https://github.com/amirfarzin/FVM_PDEsolver.jl/blob/{commit}{path}#{line}",
    sitename="FVM_PDEsolver.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
