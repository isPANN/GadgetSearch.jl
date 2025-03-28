module GadgetSearch

using JuMP, HiGHS, COPT, GLPK
using Suppressor
using Graphs, GraphIO, IterTools
using Colors, ColorSchemes
using Combinatorics
using JSON3, JSON
using Karnak, Luxor
using Random

include("dataloaders.jl")
include("utils.jl")
include("graphsearch.jl")
include("generateUDG.jl")
include("unitdiskgraph.jl")
include("visualize.jl")

export search_single_rule, search_rules
export find_maximal_independent_sets
export save_results_to_json
export read_graph_dict, read_graph
export format_truth_table
export show_rule_info, generic_rule, reconstruct_rule_id
export generate_degeneracy_cases
export load_gadget, load_unweighted_grid_gadget, load_grid_gadget
export plot_single_gadget
export extract_rule_ids, extract_graph_ids
export check_gadget

end