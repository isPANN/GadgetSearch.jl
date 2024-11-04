function loadJSONFile(filename::String)
    data = JSON.parsefile(filename)
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

function findByDegeneracy(data::Dict{Vector{String}, NamedTuple}, degeneracy::Vector{String})
    return get(data, degeneracy, "Degeneracy not found.")
end

function findByDegeneracy(filename::String, degeneracy::Vector{String})
    data = loadJSONFile(filename)
    return findByDegeneracy(data, degeneracy)
end

function findByDegeneracy(filename::String, bit_num::Int, degeneracy::Vector{Int})
    degeneracy_key = Vector{String}[join(bin(elem, bit_num), "") for elem in degeneracy]
    return findByDegeneracy(filename, degeneracy_key)
end

# function findByGateID(data::Dict{Vector{String}, NamedTuple}, gate_id::Int)
#     return get(data, gate_id, "Gate not found.")
# end

# function findByGateID(filename::String, gate_id::Int)
#     data = loadJSONFile(filename)
#     return findByGateID(data, gate_id)
# end


