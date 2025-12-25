using GadgetSearch
using Documenter
using Literate

# Literate
for each in readdir(pkgdir(GadgetSearch, "examples"))
    input_file = pkgdir(GadgetSearch, "examples", each)
    endswith(input_file, ".jl") || continue
    @info "building" input_file
    output_dir = pkgdir(GadgetSearch, "docs", "src", "generated")
    Literate.markdown(input_file, output_dir; name=each[1:end-3], execute=false)
end

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
    doctest = ("doctest=true" in ARGS),
    pages=[
        "Home" => "index.md",
        "Examples" => [
            "Rydberg Gadgets on Triangular Lattice" => "generated/trangular_Rydberg_example.md",
            "QUBO Gadgets on Triangular Lattice" => "generated/triangular_QUBO_example.md",
        ],
        "Reference" => "ref.md",
    ],
)

deploydocs(;
    repo="github.com/isPANN/GadgetSearch.jl",
    devbranch="main",
)
