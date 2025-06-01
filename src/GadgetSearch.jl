module GadgetSearch

# Core dependencies
using JuMP
using HiGHS
using Suppressor
using Graphs, GraphIO
# using HDF5
using ProgressMeter
# using JLD2
using Serialization
# using SQLite
# using Serialization
# using DBInterface
# using DataFrames
using IterTools
using Combinatorics
using LinearAlgebra
using Random
using Base.Threads

# Visualization and data handling
using Colors, ColorSchemes
using JSON3, JSON
using Karnak, Luxor

# Type definitions
# include("types/gadgets.jl")
# include("types/parameters.jl")
# include("types/graph_types.jl")

# # Include core functionality
# include("core/graphsearch.jl")
# include("core/generateudg.jl")

# # Include utilities
# include("utils/dataloaders.jl")
# include("utils/utils.jl")

# # Include settings and visualization
# include("settings.jl")
# include("visualization/visualize.jl")

# After refactoring
include("graphio/graph6.jl")
include("graphio/graphloader.jl")

include("core/search.jl")

# Export public API
# Core search functions
export search_single_rule, search_rules

# Type definitions
export SearchParameters, AbstractGadget, Gadget

# Trait system
export AbstractGraphType, GeneralGraph, GridGraph
export AbstractSearchStrategy, GenericConstraintSearch, LogicGateSearch

# Graph utilities
export find_maximal_independent_sets
export read_graph_dict, read_graph

# Data handling
export save_results_to_json
export format_truth_table
export extract_rule_ids, extract_graph_ids

# Rule utilities
export show_rule_info, generic_rule, reconstruct_rule_id
export generate_degeneracy_cases

# Gadget handling
export load_gadget, load_unweighted_grid_gadget, load_grid_gadget, load_grid_gadget_old
export check_gadget

# Visualization
export plot_single_gadget, plot_single_gadget_new, plot_single_gadget_mis

# UDG generation
export generate_grid_udgs, _split_large_file


# After refactoring
# export save_g6_graph, save_g6_graphs
# export read_g6_graph, read_g6_graphs
# export get_gids
# export each_g6_graph


export GraphDataset
export GraphLoader
export search_over_dataset


end # module