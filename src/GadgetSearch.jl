module GadgetSearch

# Core dependencies
using JuMP
# using Suppressor
using Graphs, GraphIO
using ProgressMeter
using IterTools
using Combinatorics
using Random

# Visualization and data handling
using Colors, ColorSchemes
using JSON3
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
export search_by_truth_tables
export save_cache
export save_graph

export Square, Triangular
export generate_full_grid_udg
export generate_full_grid_graph

# Core types
export Gadget
export save_results_to_json

# Energy models
export EnergyModel, RydbergModel, QUBOModel

# Constraint types
export GadgetConstraint, TruthTableConstraint, StateConstraint

# Search functions
export search_gadgets
export search_by_truth_tables
export search_by_state_constraints

# Visualization
export get_radius
export plot_gadget

# Utilities
export clear_cache!
export get_cache_stats
export check_gadget, check_gadget_rydberg, check_gadget_qubo

end # module
