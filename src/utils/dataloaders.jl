# Note: Gadget struct is now defined in GadgetSearch.jl

"""
Represents a gadget with grid position information

# Fields
- `rule_id::Int`: ID of the logic gate
- `ground_states::Vector{String}`: Ground states of the gadget
- `graph_id::Int`: ID of the graph
- `graph::SimpleGraph{Int}`: The graph structure
- `pins::Vector{Int}`: Pin vertices
- `weights::AbstractVector`: Vertex weights
- `pos::Vector{Tuple{Int, Int}}`: Grid positions of vertices
"""
struct grid_gadget{T<:Real} <: AbstractGadget
    rule_id::Int
    ground_states::Vector{String}
    graph_id::Int
    graph::SimpleGraph{Int}
    pins::Vector{Int}
    weights::Vector{T}
    pos::Vector{Tuple{Int, Int}}
end

"""
Represents an unweighted gadget with grid position information

# Fields
- `graph::SimpleGraph{Int}`: The graph structure
- `pins::Vector{Int}`: Pin vertices
- `pos::Vector{Tuple{Int, Int}}`: Grid positions of vertices
"""
struct unweighted_grid_gadget <: AbstractGadget
    graph::SimpleGraph{Int}
    pins::Vector{Int}
    pos::Vector{Tuple{Int, Int}}
end

# Get a function to process weights based on rounding preference
function _get_weight_function(round_weight::Bool)
    if round_weight
        @info "The weights will be converted to integers by rounding."
        return x -> Int.(round.(x))  # return round function
    else
        return x -> x  # return identity function
    end
end

"""
    load_gadget(gadget_file_name::String; key::Symbol=:rule_id, round_weight::Bool=true)

Load gadgets from a JSON file

# Arguments
- `gadget_file_name::String`: Path to the JSON file containing gadget data

# Keyword Arguments
- `key::Symbol=:rule_id`: The key to use for indexing the result dictionary
- `round_weight::Bool=true`: Whether to round weights to integers

# Returns
- A dictionary mapping keys to Gadget objects
"""
function load_gadget(gadget_file_name::String; key::Symbol=:rule_id, round_weight::Bool=true)
    # Parse the JSON file
    gadget_data = JSON.parsefile(gadget_file_name)
    result_dict = Dict{Int, Gadget}()
    weight_func = _get_weight_function(round_weight)

    # Process each gadget entry
    for entry in gadget_data
        # Extract key and rule ID
        dict_key = entry[String(key)]
        rule_id = Int(entry["rule_id"])

        # Extract and process node weights
        node_weights = [weight_func(node["weight"]) for node in entry["graph"]["nodes"]]

        # Create graph from edges
        g = SimpleGraph(length(node_weights))
        for edge in entry["graph"]["edges"]
            add_edge!(g, edge["source"], edge["target"])
        end

        # Extract pins and ground states
        pins = Int.(entry["pins"])
        ground_states = String.(entry["ground_states"])

        # Create and store the Gadget object
        result_dict[dict_key] = Gadget(
            rule_id,
            ground_states,
            -1,  # Default graph_id
            g,
            pins,
            node_weights
        )
    end

    return result_dict
end


function load_grid_gadget(gadget_file_name::String; key::Symbol=:rule_id, round_weight::Bool=true)
    # Parse the JSON files
    gadget_data = JSON.parsefile(gadget_file_name)
    # Initialize result dictionary and weight function
    result_dict = Dict{Int, grid_gadget}()
    weight_func = _get_weight_function(round_weight)

    # Process each gadget entry
    for entry in gadget_data
        # Extract key and rule ID
        dict_key = entry[String(key)]
        rule_id = Int(entry["rule_id"])

        # Extract and process node weights
        node_weights = [weight_func(node["weight"]) for node in entry["graph"]["nodes"]]

        # Create graph from edges
        g = SimpleGraph(length(node_weights))
        for edge in entry["graph"]["edges"]
            add_edge!(g, edge["source"], edge["target"])
        end

        # Extract pins and graph ID
        pins = Int.(entry["pins"])
        graph_id = Int(entry["graph"]["graph_id"])

        positions_data = entry["graph"]["positions"]
        pos = [(Int(p["position"][1]), Int(p["position"][2])) for p in sort(positions_data, by = p -> Int(p["id"]))]
        # Extract ground states
        ground_states = String.(entry["ground_states"])

        # Create and store the grid_gadget object
        result_dict[dict_key] = grid_gadget(
            rule_id,
            ground_states,
            graph_id,
            g,
            pins,
            node_weights,
            pos
        )
    end

    return result_dict
end


"""
    load_grid_gadget_old(m::Int, n::Int, pin_pad::Int, gadget_file_name::String, pos_file_name::String; key::Symbol=:rule_id, round_weight::Bool=true)

Load grid gadgets from JSON files

# Arguments
- `m::Int`: Number of rows in the grid
- `n::Int`: Number of columns in the grid
- `pin_pad::Int`: Padding around the grid for pins
- `gadget_file_name::String`: Path to the JSON file containing gadget data
- `pos_file_name::String`: Path to the JSON file containing position data

# Keyword Arguments
- `key::Symbol=:rule_id`: The key to use for indexing the result dictionary
- `round_weight::Bool=true`: Whether to round weights to integers

# Returns
- A dictionary mapping keys to grid_gadget objects
"""
function load_grid_gadget_old(m::Int, n::Int, pin_pad::Int, gadget_file_name::String, pos_file_name::String; key::Symbol=:rule_id, round_weight::Bool=true)
    # Parse the JSON files
    gadget_data = JSON.parsefile(gadget_file_name)
    pos_data = JSON.parsefile(pos_file_name)

    # Initialize result dictionary and weight function
    result_dict = Dict{Int, grid_gadget}()
    weight_func = _get_weight_function(round_weight)

    # Process each gadget entry
    for entry in gadget_data
        # Extract key and rule ID
        dict_key = entry[String(key)]
        rule_id = Int(entry["rule_id"])

        # Extract and process node weights
        node_weights = [weight_func(node["weight"]) for node in entry["graph"]["nodes"]]

        # Create graph from edges
        g = SimpleGraph(length(node_weights))
        for edge in entry["graph"]["edges"]
            add_edge!(g, edge["source"], edge["target"])
        end

        # Extract pins and graph ID
        pins = Int.(entry["pins"])
        graph_id = Int(entry["graph"]["graph_id"])

        # Convert position indices to tuples
        pos_int = Int.(pos_data["$graph_id"])
        pos = _index_to_tuple(pos_int, (m+2pin_pad, n+2pin_pad))

        # Extract ground states
        ground_states = String.(entry["ground_states"])

        # Create and store the grid_gadget object
        result_dict[dict_key] = grid_gadget(
            rule_id,
            ground_states,
            graph_id,
            g,
            pins,
            node_weights,
            pos
        )
    end

    return result_dict
end

"""
    load_unweighted_grid_gadget(m::Int, n::Int, pin_pad::Int, gadget_file_name::String, pos_file_name::String, pin_set::Vector{Int})
    
Load unweighted grid gadgets from files

# Arguments
- `m::Int`: Number of rows in the grid
- `n::Int`: Number of columns in the grid
- `pin_pad::Int`: Padding around the grid for pins
- `gadget_file_name::String`: Path to the file containing graph data
- `pos_file_name::String`: Path to the JSON file containing position data
- `pin_set::Vector{Int}`: Set of pins to use for all gadgets

# Returns
- A dictionary mapping graph IDs to unweighted_grid_gadget objects
"""
function load_unweighted_grid_gadget(m::Int, n::Int, pin_pad::Int, gadget_file_name::String, pos_file_name::String, pin_set::Vector{Int})
    # Load graph data and position data
    graph_data = read_graph_dict(gadget_file_name)
    pos_data = JSON.parsefile(pos_file_name)

    # Initialize result dictionary
    result_dict = Dict{Int, unweighted_grid_gadget}()

    # Process each graph
    for entry in keys(graph_data)
        # Get the graph and extract its ID
        g = graph_data[entry]
        graph_id = extract_gids(entry)

        # Convert position indices to tuples
        pos_int = Int.(pos_data["$graph_id"])
        pos = _index_to_tuple(pos_int, (m+2pin_pad, n+2pin_pad))

        # Create and store the unweighted_grid_gadget object
        result_dict[graph_id] = unweighted_grid_gadget(
            g,
            pin_set,
            pos
        )
    end

    return result_dict
end


