function readgraphdict(path::String)
    graph = loadgraphs(path, GraphIO.Graph6.Graph6Format())
    return graph
end

function readgraph(path::String, id::Int)
    graph = loadgraph(path, "graph$(id)", GraphIO.Graph6.Graph6Format())
    return graph
end

function plotgraphs(graphs, format::Symbol=:png; saved_path::String=".", saved_name::String="")
    valid_formats = [:png, :pdf]
    if format ∉ valid_formats
        @error "Unsupported format. Valid formats are: $(join(valid_formats,","))."
    end

    for (idx, gname) in enumerate(gadgets)
        filepath = joinpath("./SearchGadgets.jl/", "udg$idx.$png")
        draw(Compose.PNG(filepath, 16cm, 16cm), gplot(gname.g))
    end
end

function plotgraphs(graphs::Dict{String, Graphs.SimpleGraphs.SimpleGraph}, format::Symbol=:png; saved_path::String=".", saved_name::String="")
    valid_formats = [:png, :pdf]
    if format ∉ valid_formats
        @error "Unsupported format. Valid formats are: $(join(valid_formats,","))."
    end

    for gname in keys(graphs)
        filepath = joinpath(saved_path, "$saved_name$gname.$(String(format))")
        if format == :png
            draw(Compose.PNG(filepath, 16cm, 16cm), gplot(graphs[gname]))
        elseif format == :pdf
            draw(Compose.PDF(filepath, 16cm, 16cm), gplot(graphs[gname]))
        end
    end
end

function plotcoloredgraph(graph_info::NamedTuple, format::Symbol=:png; saved_path::String=".", name::String="graph")
    node_weights = graph_info.node_weights
    work_nodes = graph_info.work_nodes
    
    node_colors = [ colorant"blue" for _ in Graphs.vertices(graph_info.graph)]
    color_dict = Dict(
        1 => colorant"cyan",
        2 => colorant"red",
        3 => colorant"yellow",
        4 => colorant"green",
        5 => colorant"purple",
        6 => colorant"orange",
        7 => colorant"brown",
        8 => colorant"pink",
        9 => colorant"blue",
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
    valid_formats = [:png, :pdf]
    if format ∉ valid_formats
        @error "Unsupported format. Valid formats are: $(join(valid_formats,","))."
    end
    filepath = joinpath(saved_path, "$name.$(String(format))")
    if format == :png
        draw(Compose.PNG(filepath, 16cm, 16cm), p)
    elseif format == :pdf
        draw(Compose.PDF(filepath, 5cm, 3cm), p)
    end
end

function plotcoloredgraphs(graphs_info::Union{Dict{Integer, NamedTuple}, Dict{Vector{String}, NamedTuple}}, format::Symbol=:png; saved_path::String=".", name::String="graph", convert_to_latex::Bool=false)
    if !isdir(saved_path)
        @info "Creating directory: $saved_path"
        mkpath(saved_path)
    end
    for key in keys(graphs_info)
        plotcoloredgraph(graphs_info[key], format; saved_path=saved_path, name="$name$key")
    end
    if convert_to_latex
        open(joinpath(saved_path, "results.tex"), "w") do io
            write(io, """
            \\documentclass{article}
            \\usepackage{graphicx}
            \\usepackage{pdfpages}
            \\usepackage{geometry}
            \\usepackage{longtable}
            \\usepackage{array}
            \\usepackage{tikz}
            \\usepackage{xcolor}
            \\usepackage{pgfplots}
            \\pgfplotsset{compat=1.18}
            \\geometry{a4paper, margin=1in}
            \\newcolumntype{P}[1]{>{\\centering\\arraybackslash}p{#1}}

            \\begin{document}
            \\begin{center}
            Results
            \\begin{tikzpicture}[scale=0.5]
                \\draw[fill={rgb,1:red,0;green,1;blue,1}] (0, 0) rectangle (1, 1) node[midway] {1}; % cyan
                \\draw[fill={rgb,1:red,1;green,0;blue,0}] (1, 0) rectangle (2, 1) node[midway] {2}; % red
                \\draw[fill={rgb,1:red,1;green,1;blue,0}] (2, 0) rectangle (3, 1) node[midway] {3}; % yellow
                \\draw[fill={rgb,1:red,0;green,1;blue,0}] (3, 0) rectangle (4, 1) node[midway] {4}; % green
                \\draw[fill={rgb,1:red,0.5;green,0;blue,0.5}] (4, 0) rectangle (5, 1) node[midway] {5}; % purple
                \\draw[fill={rgb,1:red,1;green,0.5;blue,0}] (5, 0) rectangle (6, 1) node[midway] {6}; % orange
                % \\draw[fill={rgb,1:red,0.6;green,0.3;blue,0.1}] (6, 0) rectangle (7, 1) node[midway] {7}; % brown
                % \\draw[fill={rgb,1:red,1;green,0.75;blue,0.8}] (7, 0) rectangle (8, 1) node[midway] {8}; % pink
                % \\draw[fill={rgb,1:red,0;green,0;blue,1}] (8, 0) rectangle (9, 1) node[midway] {9}; % blue
                % \\draw[fill={rgb,1:red,1;green,0;blue,1}] (9, 0) rectangle (10, 1) node[midway] {10}; % magenta
            \\end{tikzpicture}
            \\end{center}
            
            \\begin{longtable}{|P{2cm}|P{5cm}|P{5cm}|}
            \\hline
            \\multicolumn{1}{|c|}{\\textbf{Index}} & 
            \\multicolumn{1}{c|}{\\textbf{Degeneracy}} & 
            \\multicolumn{1}{c|}{\\textbf{Graph}} \\\\ \\hline
            \\endfirsthead
            \\hline
            \\multicolumn{1}{|c|}{\\textbf{Index}} & 
            \\multicolumn{1}{c|}{\\textbf{Degeneracy}} & 
            \\multicolumn{1}{c|}{\\textbf{Graph}} \\\\ \\hline
            \\endhead
            """)

            for key in sort(collect(keys(graphs_info)))
                write(io, """
                $key & $(Vector{String}(graphs_info[key].key_info)) & \\includegraphics[width=0.2\\textwidth]{$(name)$(key).pdf} \\\\ \\hline
                """)
            end
            
            write(io, """
            \\end{longtable}
            \\end{document}
            """)
        end
    end
end

function checkgraphmis(graph_info::NamedTuple)
    g = graph_info.graph
    mis_result = find_maximal_independent_sets(g)
    mis_num = size(mis_result, 1)

    work_bits_value_vector = [[mis_result[i, j] for j in graph_info.work_nodes] for i in 1:mis_num]
    work_bits_value_string = [join(map(string, subarr)) for subarr in work_bits_value_vector]
    energy_value = [sum([graph_info.node_weights[j] * mis_result[i, j] for j in 1:nv(g)])  for i in 1:mis_num]
    max_value = maximum(energy_value)
    max_indices = findall(x -> x == max_value, energy_value)
    @info """
    All Maximal Independent States' value: $(work_bits_value_string).
    Corresponding energy values: $(energy_value).
    => Ground States for this graph: $(work_bits_value_string[max_indices]).
    """
    return work_bits_value_string, energy_value, work_bits_value_string[max_indices]
end

function genericgate(gate_id::Int, input_bits::Int, output_bits::Int)
    num_inputs = 2^input_bits
    num_outputs = 2^output_bits
    max_gateid = num_outputs^num_inputs

    if gate_id < 0 || gate_id > max_gateid
        @error("Gate ID must be between 0 and $max_gateid")
    end
    degeneracy = Int[]
    # @info "==== Gate ID: $gate_id ===="
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

function showgateinfo(gate_id::Int, input_bits::Int, output_bits::Int)
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


function generate_bitvectors(bit_num::Int, indices::Vector{Vector{Int}})::Matrix{Int}
    bit_vectors = zeros(Int, length(indices), bit_num)
    for (col, idxs) in enumerate(indices)
        for idx in idxs
            bit_vectors[col, idx] = 1  # 逐元素赋值
        end
    end
    return bit_vectors
end

format_degeneracy_input(degeneracy::Vector{Vector{Int}})::Vector{Int} = [sum(v[i] * 2^(length(v) - i) for i in 1:length(v)) for v in degeneracy]
format_degeneracy_input(degeneracy::Vector{String})::Vector{Int} = [parse(Int, s; base=2) for s in degeneracy]
format_degeneracy_output(degeneracy::Vector{Int}, bit_length::Int)::Vector{String} = [join(reverse(digits(x, base=2, pad=bit_length))) for x in degeneracy]
# format_degeneracy_output(degeneracy::Vector{Vector{Int}})::Vector{String} = [join(string.(row)) for row in degeneracy]

function get_candidate_degeneracy(index_matrix::AbstractMatrix{Int}, index::Vector{Int})::Vector{Int}
    result = index_matrix[:, index]
    num_cases = size(index_matrix, 1)
    num_index = length(index)
    decimal_values = zeros(Int, num_cases)
    for i in 1:num_cases
        decimal_val = 0
        for j in 1:num_index
            decimal_val += result[i, j] * 2^(num_index - j)
        end
        decimal_values[i] = decimal_val
    end
    return decimal_values
end
