function readGraphDictFile(path::String)
    graph = loadgraphs(path, GraphIO.Graph6.Graph6Format())
    return graph
end

function readGraphFile(path::String, id::Int)
    graph = loadgraph(path, "graph$(id)", GraphIO.Graph6.Graph6Format())
    return graph
end

function plotGraphs(graphs::Dict{String, Graphs.SimpleGraphs.SimpleGraph}, saved_path::String)
    for gname in keys(graphs)
        draw(Compose.PNG(saved_path * "$gname.png", 16cm, 16cm), gplot(graphs[gname]))
    end
end

function plotColoredGraph(graph_info::NamedTuple, saved_path::String, name::String="graph.png")
    node_weights = graph_info.node_weights
    work_nodes = graph_info.work_nodes
    
    node_colors = [ colorant"blue" for _ in Graphs.vertices(graph_info.graph)]
    color_dict = Dict(
        1 => colorant"blue",
        2 => colorant"red",
        3 => colorant"green",
        4 => colorant"yellow",
        5 => colorant"purple",
        6 => colorant"orange",
        7 => colorant"brown",
        8 => colorant"pink",
        9 => colorant"cyan",
        10 => colorant"magenta"
    )
    for (v, w) in node_weights
        node_colors[v] = color_dict[w]
    end
    node_labels = [ "" for _ in Graphs.vertices(graph_info.graph)]
    for (i, node) in enumerate(work_nodes)
        node_labels[node] = string(i)
    end
    p = gplot(graph_info.graph,
              nodefillc = node_colors,
              nodelabel = node_labels
    )
    draw(Compose.PNG(saved_path * name, 16cm, 16cm), p)
end

function checkGraphMIS(graph_info::NamedTuple)
    g = graph_info.graph
    maximalis = MaximalIS(g)
    mis_problem = GenericTensorNetwork(maximalis)
    mis_result = read_config(solve(mis_problem, ConfigsAll())[])
    mis_num = length(mis_result)

    work_bits_value_vector = [[Int(mis_result[i][j]) for j in graph_info.work_nodes] for i in 1:mis_num]
    work_bits_value_string = [join(map(string, subarr)) for subarr in work_bits_value_vector]
    energy_value = [sum([graph_info.node_weights[j] * Int(mis_result[i][j]) for j in 1:nv(g)])  for i in 1:mis_num]
    min_value = minimum(energy_value)
    min_indices = findall(x -> x == min_value, energy_value)
    @info "All Maximal Independent States' value: $(work_bits_value_string)."
    @info "Corresponding energy values: $(energy_value)."
    @info "=> Ground States for this graph: $(work_bits_value_string[min_indices])."
    return work_bits_value_string, energy_value, work_bits_value_string[min_indices]
end

function genericGate(gate_id::Int, input_bits::Int, output_bits::Int)
    num_inputs = 2^input_bits
    num_outputs = 2^output_bits
    max_gateid = num_outputs^num_inputs

    if gate_id < 0 || gate_id > max_gateid
        @error("Gate ID must be between 0 and $max_gateid")
    end
    degeneracy = Int[]
    @info "==== Gate ID: $gate_id ===="
    for input in 0:(num_inputs - 1)
        output = (gate_id >> (input * output_bits)) & (num_outputs - 1)
        input_bin_str = string(input, base=2, pad=input_bits)
        output_bin_str = string(output, base=2, pad=output_bits)
        combined_bin_str = input_bin_str * output_bin_str

        degen = parse(Int, combined_bin_str, base=2)
        push!(degeneracy, degen)
    end
    return degeneracy
end

function showGateInfo(gate_id::Int, input_bits::Int, output_bits::Int)
    num_inputs = 2^input_bits
    num_outputs = 2^output_bits
    max_gateid = num_outputs^num_inputs

    if gate_id < 0 || gate_id > max_gateid
        @error("Gate ID must be between 0 and $max_gateid")
    end
    @info "==== Gate ID: $gate_id ===="
    for input in 0:(num_inputs - 1)
        output = (gate_id >> (input * output_bits)) & (num_outputs - 1)
        println("Input: $(string(input, base=2, pad=input_bits)) -> Output: $(string(output, base=2, pad=output_bits))")
    end
end

function bin(x::Int, n::Int)::Vector{Int}
    return digits(x, base=2, pad=n) |> reverse
end

function decimal(binary_str::String)::Int
    parse(Int, binary_str; base=2)
end

function decimal(binary_array::Vector{Int})::Int
    decimal = 0
    n = length(binary_array)
    for i in 1:n
        decimal += binary_array[i] * 2^(n - i)
    end
    return decimal
end