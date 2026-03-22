using GadgetSearch
using Documenter
using Literate

# Literate
# Keep docs examples explicit so standalone pipeline scripts under `examples/`
# are not implicitly executed by docs CI.
doc_example_files = [
    "trangular_Rydberg_example.jl",
    "triangular_QUBO_example.jl",
]
output_dir = pkgdir(GadgetSearch, "docs", "src", "generated")
mkpath(output_dir)
# Remove stale generated markdown files first to avoid carrying old examples
# into docs CI.
for file in readdir(output_dir)
    endswith(file, ".md") || continue
    rm(joinpath(output_dir, file); force=true)
end
for file in doc_example_files
    input_file = pkgdir(GadgetSearch, "examples", file)
    isfile(input_file) || error("Missing docs example source: $input_file")
    @info "building" input_file
    Literate.markdown(input_file, output_dir; name=file[1:end-3], execute=false)
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
