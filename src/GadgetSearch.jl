module GadgetSearch

# Core dependencies
using JuMP, HiGHS
using Suppressor
using Graphs, GraphIO, IterTools
using Combinatorics
using LinearAlgebra, Statistics
using Random

# Visualization and data handling
using Colors, ColorSchemes
using JSON3, JSON
using Karnak, Luxor

# Type definitions

"""Abstract type for all gadget types in the system"""
abstract type AbstractGadget end

"""Represents a gadget with its properties"""
struct Gadget{T<:Real} <: AbstractGadget
    rule_id::Int
    ground_states::Vector{String}
    graph_id::Int
    graph::SimpleGraph{Int}
    pins::Vector{Int}
    weights::Vector{T}
 end

# Include other modules
include("dataloaders.jl")   # Data loading utilities
include("utils.jl")         # Utility functions
include("graphsearch.jl")   # Core search algorithms
include("generateUDG.jl")   # UDG generation
include("visualize.jl")     # Visualization tools

# Export public API
# Core search functions
export search_single_rule, search_rules

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
export load_gadget, load_unweighted_grid_gadget, load_grid_gadget
export check_gadget

# Visualization
export plot_single_gadget

# UDG generation
export generate_grid_udgs
end