# GadgetSearch

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://isPANN.github.io/GadgetSearch.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://isPANN.github.io/GadgetSearch.jl/dev/)
[![Build Status](https://github.com/isPANN/GadgetSearch.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/isPANN/GadgetSearch.jl/actions/workflows/CI.yml?query=branch%3Amain)

A Julia package for searching and analyzing computational gadgets in graph structures. This package provides tools for finding weighted graphs that can implement logic functions through their energy landscapes, particularly focusing on unit disk graphs (UDGs) and cellular automata-like computations.

## Features

- **Truth Table Search**: Search for gadgets that implement specific logic functions
- **Graph Loading and Management**: Efficient loading and caching of large graph datasets  
- **UDG Generation**: Generate unit disk graphs on square and triangular lattices
- **Visualization**: Render gadgets with weights, pins, and energy landscapes
- **Data Export**: Save results in JSON format for analysis

## Installation

```julia
using Pkg
Pkg.add("GadgetSearch")
```

## Quick Start

```julia
using GadgetSearch
using HiGHS  # For optimization

# Load a graph dataset
loader = GraphLoader("path/to/graphs.g6")

# Create truth tables to search for
truth_tables = [
    [1 0; 0 1],    # NOT gate (2 bits: input, output)
    [1 1; 0 0]     # Another logic function
]

# Search for gadgets implementing these truth tables
results, failed = search_by_truth_tables(
    loader, 
    truth_tables;
    optimizer = () -> HiGHS.Optimizer()
)

# Visualize a found gadget
if !isempty(results)
    plot_single_gadget_new(results[1], "gadget.pdf")
end
```

## Main Components

### Core Data Structures
- `Gadget`: Contains truth table, graph, pin assignments, weights, and positions
- `GraphLoader`: Efficient loader for g6-format graph files with caching
- `GraphDataset`: Container for graph codes and layout information

### Core Functions
- `search_by_truth_tables(loader, truth_tables; ...)`: Search for multiple logic functions
- `find_matching_gadget(loader; filter, limit)`: Find gadgets matching custom criteria
- `generate_full_grid_udg(lattice, nx, ny)`: Generate unit disk graphs on grids

### Graph Generation
- `Square` and `Triangular`: Lattice types for UDG generation
- `unit_disk_graph(positions, radius)`: Create UDG from positions
- `get_radius(lattice_type)`: Get default radius for lattice type

### Visualization
- `plot_single_gadget_new(gadget, path)`: Advanced gadget visualization
- `plot_graph(graph, path; pos)`: Basic graph plotting
- `check_gadget(gadget)`: Verify gadget energy landscape

### Utility Functions
- `save_results_to_json(results, path)`: Export results to JSON
- `generic_rule(rule_id, (inputs, outputs))`: Generate truth tables by ID

## Detailed Examples

### Basic Logic Gate Search

```julia
using GadgetSearch, HiGHS

# Load pre-generated graphs
loader = GraphLoader("datasets/udgs/3in1out.g6")

# Define an XOR-like function (3 inputs, 1 output)
# Truth table: output is 1 when odd number of inputs are 1
xor3_table = [
    1 0 0 0;  # 000 -> 0
    1 0 0 1;  # 001 -> 1  
    1 0 1 0;  # 010 -> 1
    1 0 1 1;  # 011 -> 0
    1 1 0 0;  # 100 -> 1
    1 1 0 1;  # 101 -> 0
    1 1 1 0;  # 110 -> 0
    1 1 1 1;  # 111 -> 1
]

# Search for gadgets implementing this function
results, failed = search_by_truth_tables(
    loader,
    [xor3_table];
    optimizer = () -> HiGHS.Optimizer(),
    max_result_num = 5,
    save_path = "xor3_results.json"
)

println("Found $(length(results)) XOR3 gadgets")
```

### Generating Custom UDG Datasets

```julia
# Generate 3x3 square lattice UDGs
generate_full_grid_udg(Square(), 3, 3; path="square_3x3.g6")

# Generate triangular lattice UDGs  
generate_full_grid_udg(Triangular(), 2, 4; path="triangular_2x4.g6")

# Load the generated graphs
loader = GraphLoader("square_3x3.g6")
println("Generated $(length(loader)) unique graphs")
```

### Custom Search with Filter Functions

```julia
# Create custom filter for specific constraints
function my_filter(graph, positions, pin_set)
    # Only accept graphs with exactly 5 vertices
    if nv(graph) != 5
        return (nothing, BitMatrix(undef, 0, 0), nothing)
    end
    
    # Use standard truth table matching
    truth_table = [1 0; 0 1]  # NOT gate
    filter_fn = make_filter(truth_table, () -> HiGHS.Optimizer(), nothing)
    return filter_fn(graph, positions, pin_set)
end

# Search with custom filter
results = find_matching_gadget(loader; filter=my_filter, limit=100)
```

### Working with Results

```julia
# Analyze a found gadget
if !isempty(results)
    gadget = results[1]
    
    # Check the energy landscape
    check_gadget(gadget)
    
    # Create visualizations
    plot_single_gadget_new(gadget, "gadget_detailed.pdf"; 
                          show_weights=true, 
                          background_grid=true)
    
    # Access gadget properties
    println("Pins: ", gadget.pins)
    println("Truth table size: ", size(gadget.ground_states))
    println("Graph has $(nv(gadget.graph)) vertices and $(ne(gadget.graph)) edges")
end
```

## Datasets

The package includes several pre-computed datasets in the `datasets/` directory:

### Logic Gates (`datasets/logic_gates/`)
- `2in2out.json`: Two input, two output logic functions
- `3in1out.json`: Three input, one output logic functions (cellular automaton rules)

### UDGs (`datasets/udgs/`)  
- `2in2out.json`: UDG candidates for 2-input, 2-output functions
- `3in1out.json`: UDG candidates for 3-input, 1-output functions

### State Constraints (`datasets/state_constraint/`)
- Various constraint configurations for different bit sizes

## Requirements

- Julia 1.10+
- An optimization solver (HiGHS.jl recommended)
- Optional: `shortg` tool from Nauty/Traces for canonical graph generation

## Advanced Usage

### Optimization Settings

```julia
# Use different optimization settings
results, failed = search_by_truth_tables(
    loader, truth_tables;
    optimizer = () -> HiGHS.Optimizer(),
    max_samples = 50,        # Limit weight enumeration samples
    allow_defect = false,    # Require exact truth table match
    connected = true         # Only search connected graphs
)
```

### Memory Management

```julia
# Enable caching for large datasets
loader = GraphLoader("large_dataset.g6"; 
                    enable_cache=true, 
                    max_cached=5000,
                    cachepath="graphs.cache")

# Save cache for future runs
save_cache(loader)
```

