using GadgetFinding
using Documenter

DocMeta.setdocmeta!(GadgetFinding, :DocTestSetup, :(using GadgetFinding); recursive=true)

makedocs(;
    modules=[GadgetFinding],
    authors="Xiwei Pan",
    sitename="GadgetFinding.jl",
    format=Documenter.HTML(;
        canonical="https://isPANN.github.io/GadgetFinding.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/isPANN/GadgetFinding.jl",
    devbranch="main",
)
