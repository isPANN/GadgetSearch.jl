module GadgetSearch

using JuMP, HiGHS
using GenericTensorNetworks
using COPT
using Suppressor
using Graphs, GraphIO, JSON, IterTools
using GraphPlot, Compose, Colors
using Combinatorics

# TODO: every exported function should have a docstring, please read the style guide before writing docstrings
export bin, decimal, plotgraphs, plotcoloredgraph, genericgate
export readgraphdict, readgraph, checkgraphmis
export search_any_constraint, search_single_constraint, check_single_constraint, search_gates, search_single_gate
export loadjsonfile, find_by_degeneracy, find_by_gateid, showgateinfo

include("utils.jl")
include("graphsearch.jl")
include("dataloaders.jl")

end
