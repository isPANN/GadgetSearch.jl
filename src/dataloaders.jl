function loadjsonfile(filename::String)
    data = JSON.parsefile(filename)
    # TODO: use data type instead of NamedTuple
    result_dict = Dict{Vector{String}, NamedTuple}()    
    for entry in data
        degeneracy_key = entry["degeneracy"]
        if haskey(entry, "gate_id")
            gate_id = entry["gate_id"]
        end
        node_weights = Dict(node["id"] => node["weight"] for node in entry["nodes"])
        g = SimpleGraph()
        for _ in keys(node_weights)
            add_vertex!(g)
        end
        for edge in entry["edges"]
            add_edge!(g, edge["source"], edge["target"])
        end
        work_nodes = entry["work_nodes"]
        # vertex_num = entry["vertex_num"]
        if haskey(entry, "gate_id")
            result_dict[degeneracy_key] = (
                graph = g,
                node_weights = node_weights,
                work_nodes = work_nodes, 
                gate_id = gate_id
            )
        else
            result_dict[degeneracy_key] = (
                graph = g,
                node_weights = node_weights,
                work_nodes = work_nodes
            )
        end
    end
    return result_dict
end

# TODO: remove degeneracy to groundstates, degeneracy is the number of states with the same energy.
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

# function find_by_gateid(data::Dict{Vector{String}, NamedTuple}, gate_id::Int)
#     return get(data, gate_id, "Gate not found.")
# end

# function find_by_gateid(filename::String, gate_id::Int)
#     data = loadJSONFile(filename)
#     return find_by_gateid(data, gate_id)
# end


