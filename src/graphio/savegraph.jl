function save_graph(g::Vector{SimpleGraph{T}}, path::String) where T
    open(path, "w") do io
        for graph in g
            println(io, GraphIO.Graph6._graphToG6String(graph))
        end
    end
    return path
end

function save_graph(g::Vector{Tuple{SimpleGraph, Vector{Tuple{Float64, Float64}}}}, path::String; g6_only::Bool=false)
    open(path, "w") do io
        for (graph, coords) in g
            g6 = GraphIO.Graph6._graphToG6String(graph)[11:end]
            if !g6_only
                coord_str = join(["($(x),$(y))" for (x, y) in coords], ";")
                println(io, g6, " ", coord_str)
            else
                println(io, g6)
            end
        end
    end
    return path
end

function save_graph(g::Vector{Tuple{AbstractString, Vector{Tuple{Float64, Float64}}}}, path::String)
    open(path, "w") do io
        for (g6, coords) in g
            println(io, g6, " ", join(["($(x),$(y))" for (x, y) in coords], ";"))
        end
    end
    return path
end