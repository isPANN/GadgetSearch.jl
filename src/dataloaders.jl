#TODO seperate input pins and output pins
struct gadget
    rule_id::Int
    ground_states::Vector{String}
    graph_id::Int  # reserved param
    graph::SimpleGraph{Int}
    pins::Vector{Int}
    weights::AbstractVector
end

struct grid_gadget
    # This struct is only used for data loading.
    rule_id::Int
    ground_states::Vector{String}
    graph_id::Int
    graph::SimpleGraph{Int}
    pins::Vector{Int}
    weights::AbstractVector
    pos::Vector{Tuple{Int, Int}}
end

struct unweighted_grid_gadget
    graph::SimpleGraph{Int}
    pins::Vector{Int}
    pos::Vector{Tuple{Int, Int}}
end

function _get_weight_function(round_weight::Bool)
    if round_weight
        @info "The weights will be converted to integers by rounding."
        return x -> Int.(round.(x))  # return round function
    else
        return x -> x  # return identity function
    end
end

function load_gadget(gadget_file_name::String; key::Symbol=:rule_id, round_weight::Bool=true)
    gadget_data = JSON.parsefile(gadget_file_name)
    result_dict = Dict{Int, gadget}()
    weight_func = _get_weight_function(round_weight)

    for entry in gadget_data
        dict_key = entry[String(key)]
        rule_id = Int.(entry["rule_id"])
        node_weights = [weight_func(node["weight"]) for node in entry["graph"]["nodes"]]
        
        g = SimpleGraph(length(node_weights))
        for edge in entry["graph"]["edges"]
            add_edge!(g, edge["source"], edge["target"])
        end

        pins = Int.(entry["pins"])
        groud_states = String.(entry["ground_states"])

        result_dict[dict_key] = (
            gadget(
                rule_id,
                groud_states,
                -1,
                g,
                pins,
                node_weights
            )
        )
    end
end


function load_grid_gadget(m::Int, n::Int, pin_pad::Int, gadget_file_name::String, pos_file_name::String; key::Symbol=:rule_id, round_weight::Bool=true)
    gadget_data = JSON.parsefile(gadget_file_name)
    pos_data = JSON.parsefile(pos_file_name)

    result_dict = Dict{Int, grid_gadget}()
    weight_func = _get_weight_function(round_weight)

    for entry in gadget_data
        dict_key = entry[String(key)]
        rule_id = Int.(entry["rule_id"])
        node_weights = [weight_func(node["weight"]) for node in entry["graph"]["nodes"]]
        
        g = SimpleGraph(length(node_weights))
        for edge in entry["graph"]["edges"]
            add_edge!(g, edge["source"], edge["target"])
        end

        pins = Int.(entry["pins"])
        graph_id = Int.(entry["graph"]["graph_id"])
        pos_int = Int.(pos_data["$graph_id"])
        pos = _index_to_tuple(pos_int, (m+2pin_pad, n+2pin_pad))
        groud_states = String.(entry["ground_states"])

        result_dict[dict_key] = (
            grid_gadget(
                rule_id,
                groud_states,
                graph_id,
                g,
                pins,
                node_weights,
                pos
            )
        )
    end
    return result_dict
end

function load_unweighted_grid_gadget(m::Int, n::Int, pin_pad::Int, gadget_file_name::String, pos_file_name::String, pin_set::Vector{Int})
    graph_data = read_graph_dict(gadget_file_name)
    pos_data = JSON.parsefile(pos_file_name)
    result_dict = Dict{Int, unweighted_grid_gadget}()

    for entry in keys(graph_data)
        g = graph_data[entry]
        graph_id = _extract_numbers(entry)
        pos_int = Int.(pos_data["$graph_id"])
        pos = _index_to_tuple(pos_int, (m+2pin_pad, n+2pin_pad))
        
        result_dict[graph_id] = (
            unweighted_grid_gadget(
                g,
                pin_set,
                pos
            )
        )
    end
    return result_dict
end

# function find_by_gateid(data::Dict{Integer, NamedTuple}, rule_id::Int)
#     return get(data, rule_id, "Gate not found.")
# end

# function find_by_gateid(filename::String, rule_id::Int)
#     data = load_gate_json(filename, :rule_id)
#     return find_by_gateid(data, rule_id)
# end


