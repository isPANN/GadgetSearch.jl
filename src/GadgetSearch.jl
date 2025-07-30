module GadgetSearch

# Core dependencies
using JuMP
# using Suppressor
using Graphs, GraphIO
using ProgressMeter
using IterTools
using Combinatorics
using LinearAlgebra
using Random

# Visualization and data handling
using Colors, ColorSchemes
using JSON3, JSON
using Karnak, Luxor


include("graphio/graph6.jl")
include("graphio/graphloader.jl")
include("graphio/savegraph.jl")
include("graphio/udg.jl")
include("utils/ruleio.jl")
include("core/search.jl")
include("utils/gadget.jl")
include("utils/visualize.jl")

export GraphDataset
export GraphLoader
export find_matching_gadget
export save_cache
export save_graph

export Square, Triangular
export generate_full_grid_udg

export Gadget, GadgetWithPos
export save_results_to_json

export get_radius

end # module