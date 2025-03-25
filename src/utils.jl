function read_graph_dict(path::String)::Dict{String, SimpleGraph{Int}}
    graph = Graphs.loadgraphs(path, GraphIO.Graph6.Graph6Format())
    return graph
end

function read_graph(path::String, id::Int)
    graph = Graphs.loadgraph(path, "graph$(id)", GraphIO.Graph6.Graph6Format())
    return graph
end

function _split_large_file(path::String, split_size::Int=700_000)::Vector{String}
    # generate a new directory under the same path to store the split files
    tmp_dir_path = _create_tmp_dir(path, split_size)
    split_file_paths = _split_large_file(path, tmp_dir_path, split_size)
    return split_file_paths
end

function _extract_numbers(s::String)
    # Extract `graph_id` from a key of a `graph_dict`, e.g. "graph1" -> 1.
    return parse.(Int, collect(eachmatch(r"\d+", s)) .|> x -> x.match)[1]
end

function _create_tmp_dir(path::String, split_size::Int)::Union{String, Vector{String}}
    # The temporary directory name is constructed as `tmp_<name>_<split_size>`, where `<name>` is the last component of the provided `path`.
    name = last(split(path, "/"))
    split_dir_path = joinpath(dirname(path), "tmp_$(name)_$(split_size)")
    if !isdir(split_dir_path)
        # If the directory already exists, the function does not modify its contents but instead lists the files inside it.
    #     file_list = readdir(split_dir_path)
    #     return map(x -> joinpath(split_dir_path, x), file_list)
    # else
        mkdir(split_dir_path)
    end
    return split_dir_path
    
end

function _split_large_file(path::String, output_path::String, num_lines::Int=700_000)::Vector{String}
    # Split a large file into smaller files with a specified number of lines. 
    # The function returns a list of paths to the split files.
    if !isdir(output_path)
        @info "Creating directory: $output_path"
        mkpath(output_path)
    end

    split_files = String[]
    file_count = 1
    output_file = joinpath(output_path, "part_$(file_count).g6")
    out_io = open(output_file, "w")
    push!(split_files, output_file)

    try
        line_count = 0
        for line in eachline(path)
            println(out_io, line)
            line_count += 1

            if line_count >= num_lines
                close(out_io)
                file_count += 1
                output_file = joinpath(output_path, "part_$(file_count).g6")
                out_io = open(output_file, "w")
                push!(split_files, output_file)
                line_count = 0
            end
        end
    finally
        close(out_io)
    end

    return split_files
end

function check_graph_mis(graph_info::NamedTuple)
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

function check_graph_mis(g::SimpleGraph{Int}, node_weights::Vector{Int}, pins::Vector{Int})
    mis_result = find_maximal_independent_sets(g)
    mis_num = size(mis_result, 1)

    work_bits_value_vector = [[mis_result[i, j] for j in pins] for i in 1:mis_num]
    work_bits_value_string = [join(map(string, subarr)) for subarr in work_bits_value_vector]
    energy_value = [sum([node_weights[j] * mis_result[i, j] for j in 1:nv(g)])  for i in 1:mis_num]
    max_value = maximum(energy_value)
    max_indices = findall(x -> x == max_value, energy_value)
    @info """
    All Maximal Independent States' value: $(work_bits_value_string).
    Corresponding energy values: $(energy_value).
    => Ground States for this graph: $(work_bits_value_string[max_indices]).
    """
    return work_bits_value_string, energy_value, work_bits_value_string[max_indices]
end


"""
    generic_gate(gate_id::Int, input_bits::Int, output_bits::Int; show_info::Bool=false)
                ::Vector{Int}

Generate the decimal vector `ground_states` of a generic `input_bits`-in-`output_bits`-out gate with the given `gate_id`.

# Arguments
- `gate_id::Int`: the ID of the gate.
- `input_bits::Int`: the number of input bits.
- `output_bits::Int`: the number of output bits.

# Keyword Arguments
- `show_info::Bool=false`: whether to show the truth table of the gate.

# Returns
The decimal vector of the `ground_states`.
"""
function generic_gate(gate_id::Int, input_bits::Int, output_bits::Int; show_info::Bool=false)::Vector{Int}
    num_inputs = 2^input_bits
    num_outputs = 2^output_bits
    max_gateid = num_outputs^num_inputs

    if gate_id < 0 || gate_id > max_gateid
        @error("Gate ID must be between 0 and $max_gateid")
    end

    ground_states = Int[]
    mask_output = num_outputs - 1
    shift_amount = output_bits

    for input in 0:(num_inputs - 1)
        output = (gate_id >> (input * shift_amount)) & mask_output
        combined_degen = (input << shift_amount) | output  # Combine input and output into one number
        push!(ground_states, combined_degen)
    end

    if show_info
        show_gate_info(gate_id, input_bits, output_bits)
    end
    return ground_states
end

"""
    show_gate_info(gate_id::Int, input_bits::Int, output_bits::Int)

Show the truth table of a generic `input_bits`-in-`output_bits`-out gate with the given `gate_id`.

# Arguments
- `gate_id::Int`: the ID of the gate.
- `input_bits::Int`: the number of input bits.
- `output_bits::Int`: the number of output bits.
"""
function show_gate_info(gate_id::Int, input_bits::Int, output_bits::Int)
    num_inputs = 2 ^ input_bits
    num_outputs = 2 ^ output_bits
    max_gateid = num_outputs ^ num_inputs

    if gate_id < 0 || gate_id > max_gateid
        @error("Gate ID must be between 0 and $max_gateid")
    end
    @info "==== Gate ID: $gate_id ===="
    
    # Create the header for the table
    header = "Input (Binary) | Output (Binary)"
    separator = "-" ^ length(header)
    println(header)
    println(separator)

    for input in 0:(num_inputs - 1)
        # Calculate output using bitwise operations
        output = (gate_id >> (input * output_bits)) & (num_outputs - 1)

        # Format input and output to fixed-width binary strings
        input_bin = string(input, base=2, pad=input_bits)
        output_bin = string(output, base=2, pad=output_bits)

        # Print the formatted result
        println("$input_bin | $output_bin")
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

"""
    format_truth_table(truth_table::AbstractMatrix)::Vector{Int}

Format the truth table into a vector of integers, row by row.
The truth table is a matrix where each row represents a binary number.

# Arguments
- `truth_table::AbstractMatrix`: a matrix where each row represents a binary number, which can be a `Matrix{Int}` or a `BitMatrix`.
"""
function format_truth_table(truth_table::AbstractMatrix)::Vector{Int}
    bit_length = size(truth_table, 2)
    max_value = 2^bit_length - 1  # Maximum value for the given bit length

    return [let value = sum(Int(row[i]) * 2^(bit_length - i) for i in 1:bit_length)
                value > max_value && throw(ArgumentError("Computed value $value exceeds the maximum representable value $max_value for bit length $bit_length."))
                value
            end for row in eachrow(truth_table)]
end

"""
    format_grstate_output(ground_states::Vector{Int}, bit_length::Int)::Vector{String}

Format the a vector of integers in to a vector of binary strings with a fixed bit length, which is easy to store in a JSON file.

# Arguments
- `ground_states::Vector{Int}`: a vector of decimal integers representing the ground states.
"""
function format_grstate_output(ground_states::Vector{Int}, bit_length::Int)::Vector{String}
    return [join(reverse(digits(x, base=2, pad=bit_length))) for x in ground_states]
end

function get_candidate_grstates(index_matrix::AbstractMatrix{Int}, index::Vector{Int})::Vector{Int}
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

function hamming_distance_matrix(index_matrix::AbstractMatrix{Int})
    k = size(index_matrix, 1)  # the number of binary numbers
    n = size(index_matrix, 2)  # the number of bits
    H = zeros(Int, k, k)

    # Compute the Hamming distance between each pair of binary numbers
    for i in 1:k
        for j in i+1:k
            H[i, j] = sum(index_matrix[i, :] .!= index_matrix[j, :]) 
            H[j, i] = H[i, j]  # the Hamming distance is symmetric
        end
    end
    return H
end

function find_most_distant(index_matrix::AbstractMatrix{Int}, top_k::Vector{Int})
    # find the most distant binary numbers
    k = size(index_matrix, 1)
    max_k = maximum(top_k)
    if max_k > k
        return collect(1:k)
    end
    H = hamming_distance_matrix(index_matrix)

    # Compute the minimum Hamming distance for each binary number
    min_distances = [minimum(H[i, [1:i-1; i+1:end]] ) for i in 1:k]
    # Sort the binary numbers by the minimum Hamming distance
    sorted_indices = sortperm(min_distances, rev=true)
    return sorted_indices[top_k]
end

function extract_gate_ids(file_path::String)
    data = JSON.parsefile(file_path)
    gate_ids = [item["gate_id"] for item in data if haskey(item, "gate_id")]
    return gate_ids
end

function extract_graph_ids(file_path::String)
    data = JSON.parsefile(file_path)
    graph_ids = [item["graph_id"] for item in data if haskey(item, "graph_id")]
    return graph_ids
end

# Convert (row, col) to a linear index
function _tuple_to_index(tups::Vector{Tuple{Int, Int}}, shape::Tuple{Int, Int})
    return [LinearIndices(shape)[CartesianIndex(tup)] for tup in tups]
end

# Convert linear index to (row, col)
function _index_to_tuple(indices::Vector{Int}, shape::Tuple{Int, Int})
    return [Tuple(CartesianIndices(shape)[idx]) for idx in indices]
end