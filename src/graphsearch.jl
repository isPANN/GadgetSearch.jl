"""
    search_single_rule(graph_path::String, bit_num::Union{Int, Vector{Int}},        ground_states::Vector{Int}; gate_id::Int = 0, max_file_size::Int = 1_000_000,  greedy::Bool=false, threshold::Int=0, max_samples::Int=0)
    ::Union{MISGadgetSolution, Nothing}

For a series of ground states, search for a single rule that can satisfy all of them. This particular function is designed for generic constraint search.

# Arguments
- `graph_path::String`: the path of the graph `*.g6` file.
- `bit_num::Union{Int, Vector{Int}}`: the number of pins, which also shows the kind of the rule to search for, e.g. `2` for a 2-bit generic constraint, `[2, 1]` for a 2-in-1-out logic gate.
- `ground_states::Vector{Int}`: the truth table of ground states to satisfy. Note: the default value of the ground states is a vector of decimal numbers, e.g. `[0, 1, 2, 3]` for binary numbers `[00, 01, 10, 11]` if we have 2-bit constraint.

# Keyword Arguments
- `gate_id::Int = 0`: the ID of the logic gate to search for. This is a reserved parameter for the logic gate search.
- `pin_set::Vector{Int} = Int[]`: the set of pins to search for. If not provided, the function will search for all possible pins.
- `max_file_size::Int = 1_000_000`: the maximum file size to process. If the file size exceeds this value, the function will split the file and process each part separately.
- `split_size::Int = 700_000`: the maximum number of rows (graphs) in each split file.
- `start_idx::Int = 0`: the starting index of the graph to search.
- `end_idx::Int = 0`: the ending index of the graph to search.
- `greedy::Bool = false`: a flag indicating whether to use a greedy approach. Defaults to `false`.
- `threshold::Int = 0`: the maximum number of MIS combinations to consider for sampling. Defaults to no sampling.
- `max_samples::Int = 0`: the maximum number of samples to generate. Defaults to no sampling.

# Returns
Returns an instance of `MISGadgetSolution` if a valid solution is found; otherwise, returns `nothing`.

# Notes
- The function ensures that the graph is connected before proceeding.
- The function iterates over all possible pin vectors and checks if they satisfy the given ground states within the maximal independent sets of the graph.
- If a valid candidate and weight are found, they are returned immediately. 
"""
function search_single_rule(graph_path::String, bit_num::Union{Int, Vector{Int}},   
                            ground_states::Vector{Int}; 
                            gate_id::Int=0, 
                            pin_set::Vector{Int} = Int[],
                            # File processing parameters
                            max_file_size_mb::Int=30, split_size::Int=700_000, 
                            # Search strategy parameters
                            start_idx::Int=0, end_idx::Int=0, 
                            greedy::Bool=false, threshold::Int=0, max_samples::Int=0
)::Union{GadgetSolution, Nothing}
    # This is a default version.
    if (filesize(graph_path) / (1024 * 1024)) > max_file_size_mb
        # Split the file into smaller parts and process each part separately.
        split_file_paths = _split_large_file(graph_path, split_size)

        for (i, file) in enumerate(split_file_paths)
            graph_dict = read_graph_dict(file)

            result = _execute_graph_search(graph_dict, bit_num, ground_states, gate_id, pin_set, i, split_size, start_idx, end_idx, greedy, threshold, max_samples)
            isnothing(result) || return result
        end
    else
        graph_dict = read_graph_dict(graph_path)

        result = _execute_graph_search(graph_dict, bit_num, ground_states, gate_id, pin_set, 0, 0, start_idx, end_idx, greedy, threshold, max_samples)
    end
    return result
end

function search_single_rule(graph_path::String, bit_num::Union{Int, Vector{Int}}, 
                            truth_table::AbstractMatrix; 
                            gate_id::Int = 0, 
                            pin_set::Vector{Int}=Int[], 
                            max_file_size_mb::Int = 30, split_size::Int=700_000, 
                            start_idx::Int=0, end_idx::Int=0, 
                            greedy::Bool=false, threshold::Int=0, max_samples::Int=0)
    # This function receives a truth table as the ground states.
    # Convert the truth table (Matrix) to a vector of decimal numbers.
    return search_single_rule(graph_path, bit_num, format_truth_table(truth_table); gate_id=gate_id, pin_set=pin_set, max_file_size_mb=max_file_size_mb, split_size=split_size, start_idx=start_idx, end_idx=end_idx, greedy=greedy, threshold=threshold, max_samples=max_samples)
end


"""
    search_single_rule(dir_path::String, bit_num::Vector{Int}, gate_id::Int; max_file_size::Int = 1_000_000)

This is a simplified version tailored for logic gate search. When the `gate_id` of the logic gate to be searched is known, this function can be called without figuring out the `ground_states`.

# Arguments
- `dir_path::String`: the path of the directory containing the graph `*.g6` files.
- `bit_num::Vector{Int}`: a length-2 vector representing the input and output bits of the logic gate.
- `gate_id::Int`: the ID of the logic gate to search for.
"""
function search_single_rule(graph_path::String, bit_num::Vector{Int}, gate_id::Int;
                            pin_set::Vector{Int}=Int[], 
                            max_file_size_mb::Int = 30, split_size::Int=700_000, 
                            start_idx::Int=0, end_idx::Int=0,
                            greedy::Bool=false, threshold::Int=0, max_samples::Int=0)
    # Only use the gate_id to search for the logic gate.
    @assert length(bit_num) == 2
    input_num = bit_num[1]; output_num = bit_num[2];
    return search_single_rule(graph_path, bit_num, generic_gate(gate_id, input_num, output_num); gate_id=gate_id, pin_set=pin_set, max_file_size_mb=max_file_size_mb, split_size=split_size, start_idx=start_idx, end_idx=end_idx, greedy=greedy, threshold=threshold, max_samples=max_samples)
end

"""
    search_single_rule(graph_dict::Dict{String, SimpleGraph{Int}}, bit_num::Union{Int, Vector{Int}}, ground_states::Vector{Int}; gate_id::Int=0, pin_set::Vector{Int}=Int[], split_idx::Int=0, split_size::Int=0, start_idx::Int=0, end_idx::Int=0, greedy::Bool=false, threshold::Int=0, max_samples::Int=0)

Searches for a single rule that satisfies given constraints for a single graph. 

This particular function is designed for generic constraint search.

# Arguments
- `graph_dict::Dict{String, SimpleGraph{Int}}`: The input graph dictionary, where the key is the graph name like `"graph100"` and the value is a `Graphs.SimpleGraph` object.
- `bit_num::Union{Int, Vector{Int}}`: The number of pins used to represent the given rule.
- `ground_states::Vector{Int}`: A vector of ground states, where each ground state is represented as a decimal value corresponding to its binary representation.

# Keyword Arguments
- `split_idx::Int=0`: The reserved index of the split file to process. Defaults to `0`.
- `split_size::Int=0`: The maximum number of rows (graphs) in each split file. Defaults to `0`. If one wants to process split files manually, please ensure that all dictionaries have the same length `split_size` except for the last one.
"""
function search_single_rule(graph_dict::Dict{String, SimpleGraph{Int}}, 
                            bit_num::Union{Int, Vector{Int}}, 
                            ground_states::Vector{Int};
                            gate_id::Int=0,
                            pin_set::Vector{Int}=Int[], 
                            split_idx::Int=0, split_size::Int=0,
                            start_idx::Int=0, end_idx::Int=0,
                            greedy::Bool=false, threshold::Int=0, max_samples::Int=0)
    # This is a default version for a `graph dict`.
    return _execute_graph_search(graph_dict, bit_num, ground_states, gate_id, pin_set, split_idx, split_size, start_idx, end_idx, greedy, threshold, max_samples)
end

function search_single_rule(graph_dict::Dict{String, SimpleGraph{Int}}, 
                            bit_num::Union{Int, Vector{Int}}, 
                            truth_table::AbstractMatrix; 
                            pin_set::Vector{Int}=Int[], 
                            split_idx::Int=0, split_size::Int=0,
                            start_idx::Int=0, end_idx::Int=0,
                            greedy::Bool=false, threshold::Int=0, max_samples::Int=0)
    # Use the truth table to search for a single rule on a `graph_dict`.
    return _execute_graph_search(graph_dict, bit_num, format_truth_table(truth_table), gate_id, pin_set, split_idx, split_size, start_idx, end_idx, greedy, threshold, max_samples)
end

function search_single_rule(graph_dict::Dict{String, SimpleGraph{Int}}, 
                            bit_num::Vector{Int}, gate_id::Int;
                            pin_set::Vector{Int}=Int[], 
                            split_idx::Int=0, split_size::Int=0,
                            start_idx::Int=0, end_idx::Int=0,
                            greedy::Bool=false, threshold::Int=0, max_samples::Int=0)
    # Only use `gate_id` on a `graph_dict`.
    @assert length(bit_num) == 2
    input_num = bit_num[1]; output_num = bit_num[2];
    return _execute_graph_search(graph_dict, bit_num, generic_gate(gate_id, input_num, output_num), gate_id, pin_set, split_idx, split_size, start_idx, end_idx, greedy, threshold, max_samples)
end

"""
    search_single_rule(graph::SimpleGraph{Int}, bit_num::Union{Int, Vector{Int}}, ground_states::Vector{Int}; pin_set::Vector{Int}=Int[], greedy::Bool=false, threshold::Int=0, max_samples::Int=0)

Searches for a single rule that satisfies given constraints for a single graph. 

This particular function is designed for generic constraint search.

# Arguments
- `graph::SimpleGraph{Int}`: The input graph, represented as a `Graphs.SimpleGraph` object.
- `bit_num::Union{Int, Vector{Int}}`: The number of pins used to represent the given rule.
- `ground_states::Vector{Int}`: A vector of ground states, where each ground state is represented as a decimal value corresponding to its binary representation.
"""
function search_single_rule(graph::SimpleGraph{Int}, bit_num::Union{Int, Vector{Int}}, 
                            ground_states::Vector{Int};
                            pin_set::Vector{Int}=Int[], 
                            greedy::Bool=false, threshold::Int=0, max_samples::Int=0)
    # This is a default version for a single graph.
    pins, weight = _search_on_single_graph(graph, bit_num, ground_states, pin_set, greedy, threshold, max_samples)
    if !isnothing(pins)
        return convert_to_result(0, graph, pins, weight, ground_state, gate_id)
    else
        @info "No valid solution found."
        return nothing
    end
end

function search_single_rule(graph::SimpleGraph{Int}, bit_num::Union{Int, Vector{Int}}, 
                            truth_table::AbstractMatrix; 
                            pin_set::Vector{Int}=Int[], 
                            greedy::Bool=false, threshold::Int=0, max_samples::Int=0)
    # Use the truth table to search for a single rule on one graph.
    return search_single_rule(graph, bit_num, format_truth_table(truth_table); pin_set=pin_set, greedy=greedy, threshold=threshold, max_samples=max_samples)
end

function search_single_rule(graph::SimpleGraph{Int}, bit_num::Vector{Int}, gate_id::Int;
                            pin_set::Vector{Int}=Int[], 
                            greedy::Bool=false, threshold::Int=0, max_samples::Int=0)
    # Only use `gate_id` on a single graph.
    @assert length(bit_num) == 2
    input_num = bit_num[1]; output_num = bit_num[2];
    return search_single_rule(graph, bit_num, generic_gate(gate_id, input_num, output_num); pin_set=pin_set, greedy=greedy, threshold=threshold, max_samples=max_samples)
end


function _execute_graph_search(graph_dict::Dict{String, SimpleGraph{Int}}, 
                                bit_num::Union{Int, Vector{Int}}, 
                                ground_state::Vector{Int}, 
                                gate_id::Int=0, 
                                pin_set::Vector{Int}=Int[], 
                                split_idx::Int=0, split_size::Int=0, 
                                start_idx::Int=0, end_idx::Int=0, 
                                greedy::Bool=false, threshold::Int=0, max_samples::Int=0)
    # Counter for the number of processed graphs
    count = split_idx == 0 ? 0 : (split_idx - 1) * split_size

    for gname in keys(graph_dict)
        count += 1
        # Skip graphs until reaching start_idx
        (start_idx > 0 && count < start_idx) && continue
        # Stop processing graphs after reaching end_idx
        (end_idx > 0 && count > end_idx) && break 
        # Print progress message every 10,000 graphs
        count % 10000 == 0 && println("Searched $count graphs...")
        
        # Perform the search on the current graph.
        pins, weight = _search_on_single_graph(graph_dict[gname], bit_num, ground_state, pin_set, greedy, threshold, max_samples)

        # If a valid solution is found, return the result.
        if !isnothing(pins)
            if split_idx == 0
                graph_id = _extract_numbers(gname)
            else
                graph_id = _extract_numbers(gname) + (split_idx - 1) * split_size
            end
            return convert_to_result(graph_id, graph_dict[gname], pins, weight, ground_state, gate_id)
        end
    end 
end

function _search_on_single_graph(graph::SimpleGraph{Int}, bit_num::Int, ground_states::Vector{Int}, pin_set::Vector{Int}, greedy::Bool, threshold::Int, max_samples::Int)
    # This is a version for generic constraints.
    # Ground states form a vector by converting binary numbers (e.g. 101) to their corresponding decimal values (e.g. 5).
    @assert length(ground_states) > 0 && maximum(ground_states) < 2^(sum(bit_num))
    # Graph must be connected.
    Graphs.is_connected(graph) || return nothing, nothing
    vertex_num = Graphs.nv(graph)

    # Find all maximal independent sets of this graph.
    mis_result, _ = find_maximal_independent_sets(graph)

    # Generate all pin vectors, i.e. the permutations of `bit_num` vertices in `pin_set`. 
    if isempty(pin_set)
        pin_set = 1:vertex_num
    end
    @assert length(pin_set) >= bit_num
    all_candidates = _generate_constraint_bit(pin_set, bit_num)
    
    # Iterate over all possible pin vectors and search the vertex weight.
    for candidate in all_candidates 
        # Check if the `ground_states` are contained in the MISs under the choice of `candidate`.
        target_mis_indices_all = _check_grstates_in_candidate(mis_result, candidate, ground_states)
        isempty(target_mis_indices_all) && continue
        
        # Iterate over all possible combinations of the target MIS indices. This is because one ground_state may be contained in different MISs.
        for combination in Iterators.product(target_mis_indices_all...)
            target_mis_set, wrong_mis_set = _split_mis_set(mis_result, vcat(combination...))
            weight = _find_weight(vertex_num, target_mis_set, wrong_mis_set)
            isempty(weight) && continue
            return candidate, weight
        end 
    end
    return nothing, nothing
end


function _search_on_single_graph(graph::SimpleGraph{Int}, bit_num::Vector{Int}, ground_states::Vector{Int}, pin_set::Vector{Int}, greedy::Bool, threshold::Int, max_samples::Int) 
    # This is a version for a logic gate.
    @assert length(bit_num) == 2
    @assert length(ground_states) > 0 && maximum(ground_states) < 2^(sum(bit_num))
    input_num = bit_num[1]; output_num = bit_num[2];

    Graphs.is_connected(graph) || return nothing, nothing
    vertex_num = Graphs.nv(graph)

    mis_result, mis_num = find_maximal_independent_sets(graph)
    mis_num < 2^(input_num) && return nothing, nothing
    if isempty(pin_set)
        # If the pin set is not provided, all vertices are considered as potential pins.
        pin_set = 1:vertex_num
        # Exclude some vertices.
        input_candidates = _generate_pin_set(graph, mis_result, input_num)
    else
        # If the pin set is provided, all selections in `pin_set` are considered valid by default.
        input_candidates = _generate_gate_input(pin_set, input_num)
    end
    isempty(input_candidates) && return nothing, nothing

    for candidate in input_candidates 
        remain_elements = setdiff(pin_set, candidate)
        for output_bits in permutations(remain_elements, output_num)
            candidate_full = vcat(candidate, output_bits)
            
            target_mis_indices_all = _check_grstates_in_candidate(mis_result, candidate_full, ground_states, greedy)
            isempty(target_mis_indices_all) && continue

            if threshold > 0 && max_samples > 0
                selected_combinations = _sample_possible_mis(target_mis_indices_all, threshold, max_samples)
            else
                selected_combinations = collect(Iterators.product(target_mis_indices_all...))
            end
            
            for combination in selected_combinations
                target_mis_set, wrong_mis_set = _split_mis_set(mis_result, vcat(combination...))
                weight = _find_weight(vertex_num, target_mis_set, wrong_mis_set)
                isempty(weight) && continue
                return candidate_full, weight
            end
        end
    end
    return nothing, nothing 
end


function _search_on_single_graph(graph::SimpleGraph{Int}, bit_num::Vector{Int}, gate_id::Int; pin_set::Vector{Int}=Int[], greedy::Bool=true, threshold::Int=0, max_samples::Int=0)
    @assert length(bit_num) == 2
    input_num = bit_num[1]; output_num = bit_num[2];
    return _search_on_single_graph(graph, bit_num, generic_gate(gate_id, input_num, output_num); pin_set=pin_set, greedy=greedy, threshold=threshold, max_samples=max_samples)
end


function _generate_pin_set(graph::SimpleGraph{Int}, mis_result::AbstractMatrix, input_num::Int)::Vector{Vector{Int}}
    # Exclude the vertices that are not contained in enough MISs.
    all_input_candidates = _generate_gate_input(_select_columns(mis_result, input_num), input_num)
    # Check if the input pins are mutually independent in the graph.
    input_candidates = Vector{Vector{Int}}()
    for subset in all_input_candidates
        is_independent = true
        for (u, v) in combinations(subset, 2)
            if has_edge(graph, u, v)
                is_independent = false
                break
            end
        end
        is_independent && push!(input_candidates, collect(subset))
    end
    return input_candidates
end

function _sample_possible_mis(full_comb::Vector{Vector{Int}}, threshold::Int, max_samples::Int)
    total_combinations = prod(length.(full_comb))
    if total_combinations > threshold
        # when the total number of combinations exceeds the threshold, randomly sample `max_samples` combinations.
        sampled_combinations = [map(rand, full_comb) for _ in 1:max_samples]
        return sampled_combinations
    else
        # otherwise, return all possible combinations.
        return collect(Iterators.product(full_comb...))
    end
end



"""
    find_maximal_independent_sets(g::SimpleGraph{Int})::Tuple{AbstractMatrix{Int}, Int}

Finds all maximal independent sets (MIS) of a given simple graph `g` and returns them in a matrix format along with the count of such sets.

# Arguments
- `g::SimpleGraph{Int}`: The input graph represented as a `SimpleGraph` object with integer vertex labels.

# Returns
- A tuple containing:
  1. `AbstractMatrix{Int}`: A matrix where each row represents a maximal independent set (MIS) in bitvector format. Each column corresponds to a vertex, and a value of `1` indicates the vertex is part of the MIS.
  2. `Int`: The total number of maximal independent sets found.

# Details
The function computes the maximal independent sets of the input graph by leveraging the equivalence between finding the maximal independent set of a graph and finding the maximal clique of its complement graph. Internally, it uses `Graphs.maximal_cliques` to compute the maximal cliques of the complement graph and then converts the result into a matrix format.
"""
function find_maximal_independent_sets(g::SimpleGraph{Int})::Tuple{AbstractMatrix{Int}, Int}
    # The problem of finding the maximal independent set of a graph is equivalent to finding the maximal clique of its complement graph.
    # `ones_vertex` is a Vector{Vector{Int}} that contains vertex indices in each MIS.
    ones_vertex = Graphs.maximal_cliques(Graphs.complement(g))
    # change the format of the result to a matrix. Each row represents a MIS. Each column represents a vertex.
    mis_result = generate_bitvectors(Graphs.nv(g), ones_vertex)
    return mis_result, length(ones_vertex)
end

function _check_grstates_in_candidate(mis_result::AbstractMatrix{Int}, candidate::Vector{Int}, grstates::Vector{Int}, greedy::Bool=false)::Vector{Vector{Int}}
    candidate_grstates = get_candidate_grstates(mis_result, candidate)
    is_subset = all(x -> x in candidate_grstates, grstates)
    is_subset || return Vector{Vector{Int}}()
    possible_choices = [findall(==(x), candidate_grstates) for x in grstates]
    if greedy 
        greedy_choice = _greedy_choose_grstates(mis_result, candidate, possible_choices; most_k=[1])
        @assert length(greedy_choice) == length(grstates)
        return greedy_choice
    else
        return possible_choices
    end
end

function _greedy_choose_grstates(mis_result::AbstractMatrix{Int}, candidate, possible_choices; most_k::Vector{Int}=[1])
    greedy_choice = []
    different_part_idx = setdiff(collect(1:size(mis_result, 2)), candidate)
    for single_rule_choices in possible_choices
        if length(single_rule_choices) <= 2
            push!(greedy_choice, single_rule_choices)
            continue
        end
        different_part = mis_result[single_rule_choices, different_part_idx]
        greedy_single_choice = find_most_distant(different_part, most_k)
        push!(greedy_choice, single_rule_choices[greedy_single_choice])
    end
    return greedy_choice
end

function _split_mis_set(mis_result::AbstractMatrix{Int}, target_indices::Vector{Int})::Tuple{AbstractMatrix{Int}, AbstractMatrix{Int}}
    mask = zeros(Bool, size(mis_result, 1))
    mask[target_indices] .= true
    target_set = mis_result[mask, :]
    wrong_set = mis_result[.!mask, :]
    return target_set, wrong_set
end

# function compute_target_diffs(target_set::AbstractMatrix)
#     n = size(target_set, 1)
#     target_diffs = Matrix{eltype(target_set)}(undef, n * (n - 1) ÷ 2, size(target_set, 2))

#     idx = 1
#     for i in 1:(n-1)
#         for j in (i+1):n
#             target_diffs[idx, :] = target_set[i, :] .- target_set[j, :]
#             idx += 1
#         end
#     end
#     return target_diffs
# end

# function compute_wrong_target_diffs(target_set::AbstractMatrix, wrong_set::AbstractMatrix)
#     wrong_target_diffs = vcat([(target_set[i, :] .- wrong_set[j, :])' for i in axes(target_set, 1), j in axes(wrong_set, 1)]...)
#     return wrong_target_diffs
# end


function _find_weight(vertex_num::Int, target_set::AbstractMatrix{Int}, wrong_set::AbstractMatrix{Int})
    model = direct_model(GLPK.Optimizer())
    set_silent(model)
    set_string_names_on_creation(model, false)

    @variable(model, x[1:vertex_num])    
    @variable(model, C)

    ϵ = 1
    @constraint(model, target_set * x .== C)  # 形成常数向量
    @constraint(model, wrong_set * x .<= C - ϵ)  # 严格不等式
    for i in 1:vertex_num
        @constraint(model, 1 <= x[i])
    end
    optimize!(model)
    
    if is_solved_and_feasible(model)
        @info "Optimization successful!"
        return [value(x[i]) for i in 1:vertex_num]
    end
    return []
end

function _select_columns(mis::AbstractMatrix{Int}, input_num::Int)::Vector{Int}
    # Compute the row-wise sum of each column in the MIS matrix.
    col_sums = sum(mis, dims=1)
    # The input vertex (pin) must be contained in enough MISs.
    valid_columns = findall(x -> x >= 2^(input_num - 1), col_sums[:])
    length(valid_columns) < input_num && return []
    return valid_columns
end

function _generate_gate_input(valid_columns::Vector{Int}, input_num::Int)::Vector{Vector{Int}}
    n = length(valid_columns)
    continuous_combs = Vector{Vector{Int}}()
    for i in 1:n
        # generate all possible continuous combinations of `input_num` vertices starting from i in `valid_columns`.
        comb = [valid_columns[(i + j - 1) % n + 1] for j in 0:input_num-1]
        push!(continuous_combs, comb)
    end
    return continuous_combs
end

function _generate_constraint_bit(valid_columns::Vector{Int}, bit_num::Int)::Vector{Vector{Int}}
    # This is a version for a generic constraint.
    return collect(permutations(valid_columns, bit_num))
end

struct Node{T}
    id::Int
    weight::T
end

struct Edge
    source::Int
    target::Int
end

struct GadgetSolution{T}
    gate_id::Int
    ground_states::Vector{String}
    graph_id::Int
    graph::SimpleGraph{Int}
    nodes::Vector{GadgetSearch.Node{T}}
    edges::Vector{GadgetSearch.Edge}
    pins::Vector{Int}
end

function convert_to_result(graph_id::Int, g::SimpleGraph{Int}, pins::Vector{Int}, weight::Vector{T}, ground_states::Vector{Int}, gate_id::Int=0) where T
    nodes = [Node(i, weight[i]) for i in 1:length(weight)]
    edges = [Edge(src(e), dst(e)) for e in Graphs.edges(g)]
    ground_state_io = format_grstate_output(ground_states, length(pins))

    @info """ === Result ===
    Gate ID: $gate_id
    Ground States: $ground_state_io
    Graph ID: $(graph_id)
    Edges: $(join([string(src(e), " -- ", dst(e)) for e in Graphs.edges(g)], ", "))
    Nodes: $(join([string(node.id, " -- ", node.weight) for node in nodes], ", "))
    """
    return GadgetSolution{T}(
        gate_id,
        ground_state_io,
        graph_id,
        g,
        nodes,
        edges,
        pins
    )
end

function save_results_to_json(results, file_path::String)
    json_data = JSON3.write(results; pretty=true)
    open(file_path, "w") do io
        write(io, json_data)
    end
    return file_path
end

