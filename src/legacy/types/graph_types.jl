"""Abstract type for graph types in the system"""
abstract type AbstractGraphType end

"""Represents a general graph type"""
struct GeneralGraph <: AbstractGraphType end

"""
Represents a grid graph type with position information.

# Fields
- `pos_data_path::Union{String, Nothing}`: Path to position data file
- `grid_dims::Tuple{Int, Int}`: Dimensions of the grid
"""
struct GridGraph <: AbstractGraphType
    pos_data_path::Union{String, Nothing}
    grid_dims::Tuple{Int, Int}

    # Constructor with path to position data
    function GridGraph(pos_data_path::Union{String, Nothing}, grid_dims::Tuple{Int, Int}=(0, 0))
        new(pos_data_path, grid_dims)
    end
end

# Default constructor with empty position data
GridGraph() = GridGraph(nothing) 