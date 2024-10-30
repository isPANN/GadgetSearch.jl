function loadJSONFile(filename::String)
    data = JSON.parsefile(filename)

    result_dict = Dict{Vector{String}, NamedTuple}()    
    for entry in data
        degeneracy_key = entry["degeneracy"]
        node_weights = Dict(node["id"] => node["weight"] for node in entry["nodes"])
        g = SimpleGraph()
        for _ in keys(node_weights)
            add_vertex!(g)
        end
        for edge in entry["edges"]
            add_edge!(g, edge["source"], edge["target"])
        end
        work_nodes = entry["work_nodes"]
        result_dict[degeneracy_key] = (
            graph = g,
            node_weights = node_weights,
            work_nodes = work_nodes
        )
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

