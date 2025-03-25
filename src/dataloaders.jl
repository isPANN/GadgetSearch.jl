struct grid_gadget{T}
    gate_id::Int
    ground_states::Vector{String}
    graph::SimpleGraph{Int}
    pins::Vector{Int}
    weights::Vector{T}
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

#TODO 拆分一下 constraint gadget 和 gate gadget
#TODO 把 JSON 里的二进制数重新处理一下
function load_grid_gadget(m::Int, n::Int, pin_pad::Int, weight_file_name::String, pos_file_name::String, key::Symbol=:gate_id; round_weight::Bool=true)
    gadget_data = JSON.parsefile(weight_file_name)
    pos_data = JSON.parsefile(pos_file_name)

    #TODO handle other type of keys
    result_dict = Dict{Int, grid_gadget}()
    weight_func = _get_weight_function(round_weight)

    for entry in gadget_data
        dict_key = entry[String(key)]
        gate_id = Int.(entry["gate_id"])
        node_weights = [weight_func(node["weight"]) for node in entry["nodes"]]
        
        g = SimpleGraph(length(node_weights))
        for edge in entry["edges"]
            add_edge!(g, edge["source"], edge["target"])
        end

        pins = Int.(entry["pins"])
        graph_id = Int.(entry["graph_id"])
        pos_int = Int.(pos_data["$graph_id"])
        pos = _index_to_tuple(pos_int, (m+2pin_pad, n+2pin_pad))
        groud_states = String.(entry["ground_states"])

        result_dict[dict_key] = (
            grid_gadget(
                gate_id,
                groud_states,
                g,
                pins,
                node_weights,
                pos
            )
        )
    end
    return result_dict
end

function load_unweighted_gadget(m::Int, n::Int, pin_pad::Int, g6_file_name::String, pos_file_name::String, pin_set::Vector{Int})
    graph_data = read_graph_dict(g6_file_name)
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


function find_by_degeneracy(data::Dict{Vector{String}, NamedTuple}, degeneracy::Vector{String})
    return get(data, degeneracy, "Degeneracy not found.")
end

function find_by_degeneracy(filename::String, degeneracy::Vector{String})
    data = load_gate_json(filename)
    return find_by_degeneracy(data, degeneracy)
end

function find_by_degeneracy(filename::String, bit_num::Int, degeneracy::Vector{Int})
    degeneracy_key = [join(bin(elem, bit_num), "") for elem in degeneracy]
    return find_by_degeneracy(filename, degeneracy_key)
end

function find_by_gateid(data::Dict{Integer, NamedTuple}, gate_id::Int)
    return get(data, gate_id, "Gate not found.")
end

function find_by_gateid(filename::String, gate_id::Int)
    data = load_gate_json(filename, :gate_id)
    return find_by_gateid(data, gate_id)
end


