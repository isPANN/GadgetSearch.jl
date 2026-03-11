module GadgetSearch

# Core dependencies
using JuMP
# using Suppressor
using Graphs, GraphIO
using ProgressMeter
using IterTools
using Combinatorics
using Random

# Tensor network dependencies for alpha tensor computation
using GenericTensorNetworks: GenericTensorNetwork, IndependentSet, SizeMax, solve
using GenericTensorNetworks: Tropical, mis_compactify!, content

# Visualization and data handling
using Colors, ColorSchemes
using JSON3
using Karnak, Luxor


include("graphio/graph6.jl")
include("graphio/graphloader.jl")
include("graphio/savegraph.jl")
include("graphio/shortg.jl")
include("graphio/lattice.jl")
include("graphio/udg.jl")
include("utils/ruleio.jl")
include("utils/flip_variants.jl")
include("core/unweighted_search.jl")
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
export generate_triangular_udg_subsets
export dedup_inner_subsets
export triangular_adjacency
export triangular_lattice_graph

# Core types
export Gadget
export save_results_to_json

# Energy models
export EnergyModel, RydbergModel, QUBOModel

# Constraint types
export GadgetConstraint, TruthTableConstraint

# Search functions
export search_gadgets
export search_by_truth_tables

# Visualization
export get_radius
export plot_gadget

# Utilities
export clear_cache!
export get_cache_stats
export check_gadget, check_gadget_rydberg, check_gadget_qubo

# Alpha tensor
export calculate_alpha_tensor
export calculate_reduced_alpha_tensor
export is_diff_by_constant
export is_gadget_replacement

# Unweighted search
export UnweightedGadget
export make_unweighted_filter
export search_unweighted_gadgets

# Multi-target unweighted search
export MultiTargetResult
export inf_mask
export pins_prefilter
export make_multi_target_filter
export search_multi_target_gadgets

# Flip variants
export generate_flip_patterns
export generate_extended_cross
export make_flip_aware_multi_target_filter
export apply_flip_to_tensor

end # module
