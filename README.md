# GadgetSearch

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://isPANN.github.io/GadgetSearch.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://isPANN.github.io/GadgetSearch.jl/dev/)
[![Build Status](https://github.com/isPANN/GadgetSearch.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/isPANN/GadgetSearch.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/isPANN/GadgetSearch.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/isPANN/GadgetSearch.jl)

# TO BE UPDATED

## [Datasets](https://github.com/isPANN/GadgetSearch.jl/tree/main/datasets)

### [Any Constraint](https://github.com/isPANN/GadgetSearch.jl/tree/main/datasets/any_constraint)

`Any constraint` means that the logic bits are unconstrained, allowing the encoding states (i.e., the ground states) to take any number and any value. For example, in the case of two logic bits under `Any constraint`, the possible ground state values can be:
```
00 # Single-fold degeneracy
01
10
11
00, 01 # Double-fold degeneracy
…
00, 01, 10 # Triple-fold degeneracy
…
…
```

This repository includes corresponding ready-to-use datasets for 2-bit and 3-bit configurations, provided in JSON file format, which can be loaded using the `loadjsonfile` function. When called, `loadjsonfile` reads the JSON file and returns a dictionary `result_dict` structured as follows:

```julia
result_dict[degeneracy_key] = (
    graph,
    node_weights,
    work_nodes
)
```

The `find_by_degeneracy` function can directly locate specific degeneracy cases and return corresponding value (i.e. a 3-tuple) in the dictionary above.

```julia
julia> using GadgetSearch

julia> dataset_path = "/datasets/any_constraint/3bits_any_constraint.json"
julia> degeneracy = ["001", "010", "101", "111"]

julia> result = find_by_degeneracy(dataset_path, degeneracy)

(graph = Graphs.SimpleGraphs.SimpleGraph{Int64}(5, [[4, 6], [5, 6], [5], [1], [2, 3], [1, 2]]), node_weights = Dict(5 => 2.0, 4 => 1.0, 6 => 1.0, 2 => 1.0, 3 => 1.0, 1 => 1.0), work_nodes = Any[2, 1, 3])
```
Once a degeneracy configuration is retrieved, `checkgraphmis` function can be used to verify all Maximal Independent States and their corresponding energy values for the graph. This function will also return all ground states associated with the given configuration.

```julia
julia> _, _, ground_states = checkgraphmis(result)

[ Info: All Maximal Independent States' value: ["010", "111", "101", "000", "001"].
[ Info: Corresponding energy values: [3.0, 3.0, 3.0, 4.0, 3.0].
[ Info: => Ground States for this graph: ["010", "111", "101", "001"].
```
For visualizing the weighted graphs, the `plotcoloredgraph` function generates a color-coded representation of the graph based on node weights.

```julia
julia> plotcoloredgraph(result, save_path)
```
### Logic Gates
In this repository, we focus on two types of logic gates:
- Two in Two out
- Three in One out

The second type is also known as a one-dimensional cellular automaton rule.

