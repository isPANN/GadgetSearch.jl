function save_graph(g::Vector{SimpleGraph{T}}, path::String) where T
    open(path, "w") do io
        for graph in g
            JSON3.write(io, Dict("g6" => graph_to_g6(graph)))
            println(io)
        end
    end
    return path
end

function save_graph(g::Vector{Tuple{SimpleGraph{T}, Vector{Tuple{Float64, Float64}}}}, path::String) where T
    open(path, "w") do io
        for (graph, coords) in g
            obj = Dict("g6" => graph_to_g6(graph), "pos" => [[x, y] for (x, y) in coords])
            JSON3.write(io, obj)
            println(io)
        end
    end
    return path
end

function save_graph(g::Vector{Tuple{S, Vector{Tuple{Float64, Float64}}}}, path::String) where S<:AbstractString
    open(path, "w") do io
        for (g6, coords) in g
            obj = Dict("g6" => String(g6), "pos" => [[x, y] for (x, y) in coords])
            JSON3.write(io, obj)
            println(io)
        end
    end
    return path
end
