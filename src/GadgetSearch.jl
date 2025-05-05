"""
GadgetSearch

A Julia package for searching and analyzing gadgets in graph structures.
Provides functionality for graph search, gadget generation, and visualization.

# Features
- Graph search algorithms
- Gadget generation and analysis
- Visualization tools
- Data loading and processing utilities
"""
module GadgetSearch

# Core dependencies
using JuMP, HiGHS
using Suppressor
using Graphs, GraphIO, IterTools
using Combinatorics
using LinearAlgebra, Statistics
using Random
using Base.Threads

# Visualization and data handling
using Colors, ColorSchemes
using JSON3, JSON
using Karnak, Luxor

# Type definitions
include("types/gadgets.jl")
include("types/parameters.jl")
include("types/graph_types.jl")

# Include core functionality
include("core/graphsearch.jl")
include("core/generateudg.jl")

# Include utilities
include("utils/dataloaders.jl")
include("utils/utils.jl")

# Include settings and visualization
include("settings.jl")
include("visualization/visualize.jl")

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
export generate_grid_udgs

end # module