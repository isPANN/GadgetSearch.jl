using GadgetSearch
using Documenter

DocMeta.setdocmeta!(GadgetSearch, :DocTestSetup, :(using GadgetSearch); recursive=true)

makedocs(;
    modules=[GadgetSearch],
    authors="Xiwei Pan",
    sitename="GadgetSearch.jl",
    format=Documenter.HTML(;
        canonical="https://isPANN.github.io/GadgetSearch.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/isPANN/GadgetSearch.jl",
    devbranch="main",
)
