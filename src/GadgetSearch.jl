module GadgetSearch

using JuMP, HiGHS
using GenericTensorNetworks
using COPT
using Suppressor
using Graphs, GraphIO, JSON, IterTools
using GraphPlot, Compose, Colors
using Combinatorics
using JSON3

include("graphio.jl")
include("utils.jl")
include("graphsearch.jl")
include("dataloaders.jl")

# Need to reduce
export bin, plotgraphs, plotcoloredgraph, plotcoloredgraphs, genericgate
export readgraphdict, readgraph, checkgraphmis
export search_any_constraint, search_single_constraint, check_singleconstraint_for_singlegraph, search_gates, search_single_gate
export loadjsonfile, find_by_degeneracy, find_by_gateid, showgateinfo
export check_single_gate, check_single_gate_traversal, check_single_gate_traversal_form, count_nonempty_lines
export g6string_to_matrix, find_maximal_cliques
export get_candidate_degeneracy
export dict_to_new_graphs,search_add_one_vertex

# After refactor
export AbstractGraph, AdjacencyMatrix, GraphsJLGraph
export generate_bitvectors
export read_g6_file
export check_singleconstraint, check_singlegate_for_singlegraph, convert_to_result, format_degeneracy_input, format_degeneracy_output, save_results_to_json
export grid_gadgets, find_maximal_independent_sets
end