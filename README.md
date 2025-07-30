# GadgetSearch

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://isPANN.github.io/GadgetSearch.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://isPANN.github.io/GadgetSearch.jl/dev/)
[![Build Status](https://github.com/isPANN/GadgetSearch.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/isPANN/GadgetSearch.jl/actions/workflows/CI.yml?query=branch%3Amain)

A Julia package for searching and analyzing gadgets in graph structures. This package provides tools for graph search, gadget generation, and visualization, with a focus on finding and analyzing gadgets that satisfy specific constraints or implement logic gates.

## Features

- **Graph Search Algorithms**: Implement both generic constraint search and logic gate search
- **Gadget Generation**: Generate and analyze gadgets with different graph types (general and grid)
- **Visualization Tools**: Visualize gadgets and their properties
- **Data Handling**: Load, save, and process graph data in various formats

## Installation

```julia
using Pkg
Pkg.add("GadgetSearch")
```

## Quick Start

```julia
using GadgetSearch

# Search for a gadget implementing a specific rule
result = search_single_rule("path/to/graph.g6", 2; truth_table=[0 1; 1 0])

# Generate grid UDGs
udgs = generate_grid_udgs(3, 3)  # Generate 3x3 grid UDGs

# Visualize a gadget
plot_single_gadget(gadget, "output.png")
```

## Main Components

### Graph Types
- `GeneralGraph`: For general graph structures
- `GridGraph`: For grid-based graph structures with position information

### Search Strategies
- `GenericConstraintSearch`: For searching gadgets with general constraints
- `LogicGateSearch`: For searching gadgets that implement logic gates

### Core Functions
- `search_single_rule`: Search for a gadget implementing a specific rule
- `search_rules`: Search for multiple rules
- `find_maximal_independent_sets`: Find all maximal independent sets in a graph
- `generate_grid_udgs`: Generate grid-based UDGs

### Data Handling
- `load_gadget`: Load a gadget from file
- `save_results_to_json`: Save search results to JSON
- `format_truth_table`: Format truth table data

### Visualization
- `plot_single_gadget`: Visualize a single gadget

## Examples

### Searching for a Logic Gate
```julia
# Search for a 2-input, 1-output logic gate
params = SearchParameters(
    pin_set = [1, 2, 3],
    max_file_size_mb = 1,
    split_size = 1000
)

result = search_single_rule(
    "path/to/graph.g6",
    [2, 1],  # 2 inputs, 1 output
    ground_states = [0, 3],  # XOR-like behavior
    params = params
)
```

### Working with Grid Graphs
```julia
# Create a grid graph with position data
grid = GridGraph("positions.json", (3, 3))

# Search on the grid graph
result = search_single_rule(
    "path/to/graph.g6",
    2,
    graph_type = grid
)
```


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

