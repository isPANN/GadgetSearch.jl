function readGraphDictFile(path::String)
    graph = loadgraphs(path, GraphIO.Graph6.Graph6Format())
    return graph
end

function readGraphFile(path::String, id::Int)
    graph = loadgraph(path, "graph$(id)", GraphIO.Graph6.Graph6Format())
    return graph
end

function save_graph_to_json(g::AbstractGraph, filename::String)
    nodes = [Dict("id" => v) for v in vertices(g)]
    edges = [Dict("source" => src(e), "target" => dst(e)) for e in edges(g)]
    
    graph_data = Dict("nodes" => nodes, "edges" => edges)
    open(filename, "w") do file
        write(file, JSON.json(graph_data))
    end
    println("Graph saved to $filename")
end

function plotGraph(graphs::Dict{String, Graphs.SimpleGraphs.SimpleGraph}, saved_path::String="../data/graphplot/")
    for gname in keys(graphs)
        draw(Compose.PNG(saved_path * "$gname.png", 16cm, 16cm), gplot(graphs[gname]))
    end
end

# function generic2bitGate(input1::Bool, input2::Bool, gate::Symbol)
#     operations = Dict(
#         :AND => (x, y) -> x & y,
#         :OR  => (x, y) -> x | y,
#         :NOR => (x, y) -> !(x | y),
#         :XOR => (x, y) -> x ⊻ y
#     )
#     return operations[gate](input1, input2)
# end

# function generic2bitGate(input1::Int, input2::Int, ruleid::Int)
#     return (ruleid >> (input1 << 1 | input2)) & 1
# end

# function genericCAGate(left::Int, middle::Int, right::Int, ruleid::Int)
#     # N: CA rule id, 0 <= N <= 255
#     # p, q, r: input states, 0 or 1
#     return (ruleid >> (left << 2 | middle << 1 | right)) & 1
# end


function bin(x::Int, n::Int)
    return digits(x, base=2, pad=n) |> reverse
end

function decimal(binary_array::Vector{Int})
    decimal = 0
    n = length(binary_array)
    for i in 1:n
        decimal += binary_array[i] * 2^(n - i)
    end
    return decimal
end