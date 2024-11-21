function loadjsonfile(filename::String, key::Symbol=:degeneracy)
    valid_keys = [:degeneracy, :gate_id]
    if !(key in valid_keys)
        error("Invalid key. Valid keys are: $(join(valid_keys, ","))")
    end

    the_other_key = key == :degeneracy ? :gate_id : :degeneracy

    data = JSON.parsefile(filename)
    result_dict = key == :degeneracy ? Dict{Vector{String}, NamedTuple}() : Dict{Integer, NamedTuple}()

    for entry in data
        dict_key = entry[String(key)]
        node_weights = Dict(node["id"] => node["weight"] for node in entry["nodes"])

        g = SimpleGraph(length(node_weights))
        for edge in entry["edges"]
            add_edge!(g, edge["source"], edge["target"])
        end

        work_nodes = entry["work_nodes"]

        key_info = get(entry, String(the_other_key), nothing)
        result_dict[dict_key] = (
            graph = g,
            node_weights = node_weights,
            work_nodes = work_nodes,
            key_info = key_info
        )
    end
    return result_dict
end

function find_by_degeneracy(data::Dict{Vector{String}, NamedTuple}, degeneracy::Vector{String})
    return get(data, degeneracy, "Degeneracy not found.")
end

function find_by_degeneracy(filename::String, degeneracy::Vector{String})
    data = loadjsonfile(filename)
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
    data = loadjsonfile(filename, :gate_id)
    return find_by_gateid(data, gate_id)
end


