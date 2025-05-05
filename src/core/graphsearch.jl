# Search parameters are defined in types.jl

"""_process_ground_states(ground_states::Vector{Int}, truth_table::AbstractMatrix, bit_num::Union{Int, Vector{Int}}, rule_id::Int)

Process ground states from either direct input, truth table, or rule ID.

# Arguments
- `ground_states::Vector{Int}`: Ground states to satisfy
- `truth_table::AbstractMatrix`: Truth table of the logic gate
- `bit_num::Union{Int, Vector{Int}}`: Number of bits or vector of input/output bits
- `rule_id::Int`: ID of the logic gate

# Returns
- Tuple of (ground_states, rule_id)
"""
function _process_ground_states(ground_states::Vector{Int}, truth_table::AbstractMatrix, bit_num::Union{Int, Vector{Int}}, rule_id::Int)
    if isempty(ground_states)
        if isempty(truth_table)
            if length(bit_num) == 2 && rule_id >= 0
                # When the `rule_id` of the logic gate to be searched is known, this function can be called without figuring out the `ground_states`.
                ground_states = generic_rule(rule_id, bit_num)
            elseif length(bit_num) == 1 && rule_id > 0
                # When the `rule_id` of the state constraint to be searched is known, this function can be called without figuring out the `ground_states`.
                ground_states = generic_rule(rule_id, bit_num)
            else
                @error "The rule id $rule_id is not valid."
                return nothing, nothing
            end
        else
            # This function receives a truth table as the ground states.
            ground_states = format_truth_table(truth_table)
        end
    end
    if length(bit_num) == 2 && rule_id < 0
        rule_id = reconstruct_rule_id(ground_states, bit_num)
    end
    return ground_states, rule_id
end


"""
    search_rules(graph_path, bit_num, gate_list, save_path;
                    pin_set::Vector{Int} = Int[],
                    # File processing parameters
                    max_file_size_mb::Int=30, split_size::Int=700_000,
                    # Search strategy parameters
                    start_idx::Int=0, end_idx::Int=0,
                    greedy::Bool=false, threshold::Int=0, max_samples::Int=0,
                    # Grid graph parameters
                    is_grid_graph::Bool=false, pos_data_path=nothing, grid_dims=(0,0),
                    optimizer=HiGHS.Optimizer, env=nothing)

# Arguments
- `graph_path`: Path to the graph database
- `bit_num`: Number of bits or vector of input/output bits
- `gate_list`: List of gate IDs to search for
- `save_path`: Path to save results

# Keyword Arguments
- `optimizer=HiGHS.Optimizer`: Optimization solver to use
- Other parameters passed to SearchParameters constructor

# Returns
- List of rule IDs for which no solution was found
"""
function search_rules(graph_path, bit_num, gate_list, save_path;
                    pin_set::Vector{Int} = Int[],
                    # File processing parameters
                    max_file_size_mb::Int=30, split_size::Int=700_000,
                    # Search strategy parameters
                    start_idx::Int=0, end_idx::Int=0,
                    greedy::Bool=false, threshold::Int=0, max_samples::Int=0,
                    # Grid graph parameters
                    is_grid_graph::Bool=false, pos_data_path=nothing, grid_dims=(0,0),
                    # Optimizer Setttings
                    optimizer=HiGHS.Optimizer, env=nothing)

    # Determine graph type using traits
    graph_type = is_grid_graph ? GridGraph(pos_data_path, grid_dims) : GeneralGraph()

    # Create appropriate result vector based on graph type
    save_res = create_result_vector(graph_type)

    # Create search parameters object
    params = SearchParameters(
        pin_set=pin_set,
        max_file_size_mb=max_file_size_mb,
        split_size=split_size,
        start_idx=start_idx,
        end_idx=end_idx,
        greedy=greedy,
        threshold=threshold,
        max_samples=max_samples
    )

    # Prepare graph chunks for processing
    graph_chunks = if (filesize(graph_path) / (1024 * 1024)) > params.max_file_size_mb
        # Split the file into smaller parts and return as (index, path) pairs
        [(i, path) for (i, path) in enumerate(_split_large_file(graph_path, params.split_size))]
    else
        # Just use the original file with index 0
        [(0, graph_path)]
    end

    # Track rules for which solutions have not yet been found
    remaining_rules = Set(gate_list)

    # Determine search strategy
    search_strategy = determine_search_strategy(bit_num)

    # Pre-compute ground states for all rules
    rule_ground_states = Dict(rule_id => generic_rule(rule_id, bit_num) for rule_id in gate_list)

    # Process each graph chunk
    for (chunk_idx, chunk_path) in graph_chunks
        # Exit early if solutions for all rules have been found
        isempty(remaining_rules) && break

        # Load graph dictionary once per chunk
        graph_dict = read_graph_dict(chunk_path)

        # Process each graph
        for (i, (gname, graph)) in enumerate(graph_dict)
            # Skip if graph is not connected
            Graphs.is_connected(graph) || continue

            # Calculate original graph ID
            graph_id = _extract_numbers(gname) + (chunk_idx == 0 ? 0 : (chunk_idx - 1) * params.split_size)

            # Find maximal independent sets (MIS) - this is a computationally intensive operation, only needs to be done once
            mis_result, mis_num = find_maximal_independent_sets(graph)

            # Print progress message
            i % 1000 == 0 && println("Searched $(i)-th graph...")

            # Choose different processing methods based on search strategy type
            if isa(search_strategy, LogicGateSearch)
                input_num = bit_num[1]
                output_num = bit_num[2]

                # Skip current graph if there aren't enough MISs
                if mis_num < 2^(input_num)
                    continue
                end

                # Generate input pin candidates
                local pin_set = params.pin_set
                if isempty(pin_set)
                    # If no pin set is provided, all vertices are considered potential pins
                    pin_set = collect(1:Graphs.nv(graph))
                    # Generate valid input pin combinations
                    input_candidates = _generate_pin_set(graph, mis_result, input_num)
                else
                    # If pin set is provided, all choices in pin_set are considered valid by default
                    input_candidates = _generate_gate_input(pin_set, input_num)
                end

                # Skip current graph if there are no valid input candidates
                isempty(input_candidates) && continue

                # Try each input candidate
                for candidate in input_candidates
                    # Find remaining elements for output pins
                    remain_elements = setdiff(pin_set, candidate)

                    # Try each possible output pin combination
                    for output_bits in permutations(remain_elements, output_num)
                        # Combine input and output pins
                        candidate_full = vcat(candidate, output_bits)

                        # Process each remaining rule
                        for rule_id in collect(remaining_rules)
                            # Get ground states for this rule
                            ground_states = rule_ground_states[rule_id]

                            # Check if ground states are contained in MIS
                            target_mis_indices_all = _check_grstates_in_candidate(mis_result, candidate_full, ground_states, params.greedy)
                            isempty(target_mis_indices_all) && continue

                            # Find weights
                            weight = _find_weight_new(mis_result, target_mis_indices_all, optimizer, env)
                            isempty(weight) && continue

                            # If solution is found
                            result = convert_to_gadget(graph_type, graph_id, graph, candidate_full, weight, ground_states, rule_id)

                            # Save result
                            push!(save_res, result)
                            save_results_to_json(save_res, save_path)

                            # Remove from remaining rules
                            delete!(remaining_rules, rule_id)

                            # Exit early if solutions for all rules have been found
                            isempty(remaining_rules) && break
                        end

                        # Exit early if solutions for all rules have been found
                        isempty(remaining_rules) && break
                    end

                    # Exit early if solutions for all rules have been found
                    isempty(remaining_rules) && break
                end
            else  # For GenericConstraintSearch strategy
                # Generate all pin vectors, i.e., permutations of `bit_num` vertices from `pin_set`
                local pin_set = params.pin_set
                if isempty(pin_set)
                    pin_set = collect(1:Graphs.nv(graph))
                end
                @assert length(pin_set) >= bit_num
                all_candidates = _generate_constraint_bit(pin_set, bit_num)

                # Iterate through all possible pin vectors
                for candidate in all_candidates
                    # Process each remaining rule
                    for rule_id in collect(remaining_rules)
                        # Get ground states for this rule
                        ground_states = rule_ground_states[rule_id]

                        # Check if `ground_states` are contained in MIS when `candidate` is selected
                        target_mis_indices_all = _check_grstates_in_candidate(mis_result, candidate, ground_states, params.greedy)
                        isempty(target_mis_indices_all) && continue

                        # Find weights
                        weight = _find_weight_new(mis_result, target_mis_indices_all, optimizer, env)
                        isempty(weight) && continue

                        # If solution is found
                        result = convert_to_gadget(graph_type, graph_id, graph, candidate, weight, ground_states, rule_id)

                        # Save result
                        push!(save_res, result)
                        save_results_to_json(save_res, save_path)

                        # Remove from remaining rules
                        delete!(remaining_rules, rule_id)

                        # Exit early if solutions for all rules have been found
                        isempty(remaining_rules) && break
                    end

                    # Exit early if solutions for all rules have been found
                    isempty(remaining_rules) && break
                end
            end

            # Exit early if solutions for all rules have been found
            isempty(remaining_rules) && break
        end
    end

    # Return rules for which no solution was found
    return collect(remaining_rules)
end


"""
    search_single_rule(graph_path::String, bit_num::Union{Int, Vector{Int}};
                            ground_states::Vector{Int}=Int[],
                            truth_table::AbstractMatrix=Matrix{Int}(undef, 0, 0),
                            rule_id::Int=-1,
                            pin_set::Vector{Int} = Int[],
                            # File processing parameters
                            max_file_size_mb::Int=30, split_size::Int=700_000,
                            # Search strategy parameters
                            start_idx::Int=0, end_idx::Int=0,
                            greedy::Bool=false, threshold::Int=0, max_samples::Int=0,
                            # Graph type parameter
                            graph_type::AbstractGraphType=GeneralGraph(),
                            optimizer=HiGHS.Optimizer, env=nothing)::Union{AbstractGadget, Nothing}

# Arguments
- `graph_path::String`: Path to the graph database
- `bit_num::Union{Int, Vector{Int}}`: Number of bits or vector of input/output bits

# Keyword Arguments
- `ground_states::Vector{Int}=Int[]`: Ground states to satisfy
- `truth_table::AbstractMatrix=Matrix{Int}(undef, 0, 0)`: Truth table of the logic gate
- `rule_id::Int=-1`: ID of the logic gate
- `optimizer=HiGHS.Optimizer`: Optimization solver to use
- Other parameters passed to SearchParameters constructor

# Returns
- A Gadget object if a solution is found, otherwise nothing
"""
function search_single_rule(graph_path::String, bit_num::Union{Int, Vector{Int}};
                            ground_states::Vector{Int}=Int[],
                            truth_table::AbstractMatrix=Matrix{Int}(undef, 0, 0),
                            rule_id::Int=-1,
                            pin_set::Vector{Int} = Int[],
                            # File processing parameters
                            max_file_size_mb::Int=30, split_size::Int=700_000,
                            # Search strategy parameters
                            start_idx::Int=0, end_idx::Int=0,
                            greedy::Bool=false, threshold::Int=0, max_samples::Int=0,
                            # Graph type parameter
                            graph_type::AbstractGraphType=GeneralGraph(),
                            optimizer=HiGHS.Optimizer, env=nothing
)::Union{AbstractGadget, Nothing}
    # Create search parameters object
    params = SearchParameters(
        pin_set=pin_set,
        max_file_size_mb=max_file_size_mb,
        split_size=split_size,
        start_idx=start_idx,
        end_idx=end_idx,
        greedy=greedy,
        threshold=threshold,
        max_samples=max_samples
    )

    # Process ground states
    ground_states, rule_id = _process_ground_states(ground_states, truth_table, bit_num, rule_id)

    # Prepare graph chunks for processing
    graph_chunks = if (filesize(graph_path) / (1024 * 1024)) > params.max_file_size_mb
        # Split the file into smaller parts and return as (index, path) pairs
        [(i, path) for (i, path) in enumerate(_split_large_file(graph_path, params.split_size))]
    else
        # Just use the original file with index 0
        [(0, graph_path)]
    end

    # Process each graph chunk until a solution is found
    for (chunk_idx, chunk_path) in graph_chunks
        # Load graph dictionary once per chunk
        graph_dict = read_graph_dict(chunk_path)

        # Execute search on this chunk of the graph database
        result = _execute_graph_search(
            graph_dict, bit_num, ground_states, params, chunk_idx, rule_id,
            optimizer, env;
            graph_type=graph_type
        )

        # Return early if a solution is found
        isnothing(result) || return result
    end

    # No solution found
    return nothing
end


"""
    search_single_rule(graph_dict::Dict{String, SimpleGraph{Int}}, bit_num::Union{Int, Vector{Int}}; ground_states::Vector{Int}, truth_table::AbstractMatrix=Matrix{Int}(undef, 0, 0), rule_id::Int=0, pin_set::Vector{Int}=Int[], split_idx::Int=0, split_size::Int=0, start_idx::Int=0, end_idx::Int=0, greedy::Bool=false, threshold::Int=0, max_samples::Int=0)

Searches for a single rule that satisfies given constraints for a single graph.

This particular function is designed for generic constraint search.

# Arguments
- `graph_dict::Dict{String, SimpleGraph{Int}}`: Dictionary of graphs
- `bit_num::Union{Int, Vector{Int}}`: Number of bits or vector of input/output bits

# Keyword Arguments
- `ground_states::Vector{Int}=Int[]`: Ground states to satisfy
- `truth_table::AbstractMatrix=Matrix{Int}(undef, 0, 0)`: Truth table of the logic gate
- `rule_id::Int=-1`: ID of the logic gate
- `optimizer=HiGHS.Optimizer`: Optimization solver to use
- `split_idx::Int=0`: Index of the split file
- `split_size::Int=0`: Size of each split file
- Other parameters passed to SearchParameters constructor

# Returns
- A Gadget object if a solution is found, otherwise nothing
"""
function search_single_rule(graph_dict::Dict{String, SimpleGraph{Int}},
                            bit_num::Union{Int, Vector{Int}};
                            ground_states::Vector{Int}=Int[],
                            truth_table::AbstractMatrix=Matrix{Int}(undef, 0, 0),
                            rule_id::Int=-1,
                            pin_set::Vector{Int}=Int[],
                            split_idx::Int=0, split_size::Int=700_000,
                            start_idx::Int=0, end_idx::Int=0,
                            greedy::Bool=false, threshold::Int=0, max_samples::Int=0,
                            graph_type::AbstractGraphType=GeneralGraph(),
                            optimizer=HiGHS.Optimizer, env=nothing)
    # Create search parameters object
    params = SearchParameters(
        pin_set=pin_set,
        max_file_size_mb=30, # Not used for this version but needed for struct
        split_size=split_size,
        start_idx=start_idx,
        end_idx=end_idx,
        greedy=greedy,
        threshold=threshold,
        max_samples=max_samples
    )

    # Process ground states
    ground_states, rule_id = _process_ground_states(ground_states, truth_table, bit_num, rule_id)

    # Execute search on the graph dictionary
    return _execute_graph_search(
        graph_dict, bit_num, ground_states, params, split_idx, rule_id,
        optimizer, env;
        graph_type=graph_type
    )
end

"""
    search_single_rule(graph::SimpleGraph{Int}, bit_num::Union{Int, Vector{Int}}; ground_states::Vector{Int}, truth_table::AbstractMatrix=Matrix{Int}(undef, 0, 0), rule_id::Int=0, pin_set::Vector{Int}=Int[], greedy::Bool=false, threshold::Int=0, max_samples::Int=0)

Searches for a single rule that satisfies given constraints for a single graph.

This particular function is designed for generic constraint search.

# Arguments
- `graph::SimpleGraph{Int}`: The graph to search
- `bit_num::Union{Int, Vector{Int}}`: Number of bits or vector of input/output bits

# Keyword Arguments
- `ground_states::Vector{Int}=Int[]`: Ground states to satisfy
- `truth_table::AbstractMatrix=Matrix{Int}(undef, 0, 0)`: Truth table of the logic gate
- `rule_id::Int=-1`: ID of the logic gate
- `optimizer=HiGHS.Optimizer`: Optimization solver to use
- Other parameters passed to SearchParameters constructor

# Returns
- A Gadget object if a solution is found, otherwise nothing
"""
function search_single_rule(graph::SimpleGraph{Int}, bit_num::Union{Int, Vector{Int}};
                            ground_states::Vector{Int}=Int[],
                            truth_table::AbstractMatrix=Matrix{Int}(undef, 0, 0),
                            rule_id::Int=-1,
                            pin_set::Vector{Int}=Int[],
                            greedy::Bool=false, threshold::Int=0, max_samples::Int=0,
                            graph_type::AbstractGraphType=GeneralGraph(),
                            optimizer=HiGHS.Optimizer, env=nothing)
    # Create search parameters object
    params = SearchParameters(
        pin_set=pin_set,
        max_file_size_mb=30, # Not used for this version but needed for struct
        split_size=700_000,  # Not used for this version but needed for struct
        start_idx=0,         # Not used for this version but needed for struct
        end_idx=0,           # Not used for this version but needed for struct
        greedy=greedy,
        threshold=threshold,
        max_samples=max_samples
    )

    # Process ground states
    ground_states, rule_id = _process_ground_states(ground_states, truth_table, bit_num, rule_id)

    # Determine search strategy based on bit_num
    search_strategy = determine_search_strategy(bit_num)

    # Search on the single graph
    pins, weight = search_on_single_graph(search_strategy, graph, bit_num, ground_states, params, optimizer, env)

    # Return result if found
    if !isnothing(pins)
        return convert_to_gadget(graph_type, 0, graph, pins, weight, ground_states, rule_id)
    else
        @info "No valid solution found in this graph."
        return nothing
    end
end


function _execute_graph_search(graph_dict::Dict{String, SimpleGraph{Int}},
                                bit_num::Union{Int, Vector{Int}},
                                ground_state::Vector{Int},
                                params::SearchParameters,
                                split_idx::Int=0,
                                rule_id::Int=-1,
                                optimizer=HiGHS.Optimizer, env=nothing;
                                graph_type::AbstractGraphType=GeneralGraph())
    # Calculate base offset for graph IDs
    base_offset = split_idx == 0 ? 0 : (split_idx - 1) * params.split_size

    # Filter and enumerate the graphs to process
    filtered_graphs = Iterators.enumerate(graph_dict) |>
        iter -> Iterators.filter(((i, _),) ->
            (params.start_idx <= 0 || base_offset + i >= params.start_idx) &&
            (params.end_idx <= 0 || base_offset + i <= params.end_idx), iter)

    # Process each graph
    for (i, (gname, graph)) in filtered_graphs
        # Print progress message every 10,000 graphs
        # (base_offset + i) % 10000 == 0 && println("Searched $(base_offset + i) graphs...")

        # Determine search strategy based on bit_num
        search_strategy = determine_search_strategy(bit_num)

        # Perform the search on the current graph
        pins, weight = search_on_single_graph(search_strategy, graph, bit_num, ground_state, params, optimizer, env)

        # If a valid solution is found, return the result
        if !isnothing(pins)
            # Calculate graph ID based on split information
            graph_id = _extract_numbers(gname) + (split_idx == 0 ? 0 : base_offset)
            return convert_to_gadget(graph_type, graph_id, graph, pins, weight, ground_state, rule_id)
        end
    end

    # No solution found
    return nothing
end


"""
    _search_on_single_graph(graph::SimpleGraph{Int}, bit_num::Vector{Int}, rule_id::Int; pin_set::Vector{Int}=Int[], greedy::Bool=true, threshold::Int=0, max_samples::Int=0, optimizer=HiGHS.Optimizer)

Search for a solution on a single graph for a specific rule ID.

# Arguments
- `graph::SimpleGraph{Int}`: The graph to search
- `bit_num::Vector{Int}`: Vector of [input_bits, output_bits]
- `rule_id::Int`: ID of the logic gate

# Keyword Arguments
- `pin_set::Vector{Int}=Int[]`: Set of pins to search for
- `greedy::Bool=true`: Whether to use greedy approach
- `threshold::Int=0`: Maximum number of MIS combinations to consider
- `max_samples::Int=0`: Maximum number of samples to generate
- `optimizer`: Optimization solver to use

# Returns
- Tuple of (pins, weights) if a solution is found, otherwise (nothing, nothing)
"""
function _search_on_single_graph(graph::SimpleGraph{Int}, bit_num::Vector{Int}, rule_id::Int; pin_set::Vector{Int}=Int[], greedy::Bool=true, threshold::Int=0, max_samples::Int=0, optimizer=HiGHS.Optimizer, env=nothing)
    @assert length(bit_num) == 2

    # Create a minimal SearchParameters object with only the needed parameters
    params = SearchParameters(
        pin_set=pin_set,
        greedy=greedy,
        threshold=threshold,
        max_samples=max_samples,
        # Default values for unused parameters
        max_file_size_mb=30,
        split_size=700_000,
        start_idx=0,
        end_idx=0
    )

    # Get ground states for the rule
    ground_states = generic_rule(rule_id, bit_num)

    # Determine search strategy and execute search
    search_strategy = determine_search_strategy(bit_num)
    return search_on_single_graph(search_strategy, graph, bit_num, ground_states, params, optimizer, env)
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

"""
    _check_grstates_in_candidate(mis_result::AbstractMatrix{Int}, candidate::Vector{Int}, grstates::Vector{Int}, greedy::Bool=false)::Vector{Vector{Int}}

Check if ground states are contained in the candidate MISs.

# Arguments
- `mis_result::AbstractMatrix{Int}`: Matrix of MIS bit vectors
- `candidate::Vector{Int}`: Candidate pin vertices
- `grstates::Vector{Int}`: Ground states to check
- `greedy::Bool=false`: Whether to use greedy approach

# Returns
- Vector of vectors of indices of MISs that contain each ground state
"""
function _check_grstates_in_candidate(mis_result::AbstractMatrix{Int}, candidate::Vector{Int}, grstates::Vector{Int}, greedy::Bool=false)::Vector{Vector{Int}}
    # Get the decimal values of the candidate pins in each MIS
    candidate_grstates = get_candidate_grstates(mis_result, candidate)

    # Check if all ground states are contained in the candidate MISs
    is_subset = all(x -> x in candidate_grstates, grstates)
    is_subset || return Vector{Vector{Int}}()

    # Find all MISs that contain each ground state
    possible_choices = [findall(==(x), candidate_grstates) for x in grstates]

    # Use greedy approach if requested
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

"""
    _find_weight(vertex_num::Int, target_set::AbstractMatrix{Int}, wrong_set::AbstractMatrix{Int}, optimizer=HiGHS.Optimizer)

Find vertex weights that satisfy the given constraints.

# Arguments
- `vertex_num::Int`: Number of vertices in the graph
- `target_set::AbstractMatrix{Int}`: Matrix of target MIS bit vectors
- `wrong_set::AbstractMatrix{Int}`: Matrix of wrong MIS bit vectors
- `optimizer`: Optimization solver to use (default: HiGHS.Optimizer)

# Returns
- Vector of vertex weights if a solution is found, otherwise empty vector
"""
function _find_weight(vertex_num::Int, target_set::AbstractMatrix{Int}, wrong_set::AbstractMatrix{Int}, optimizer=HiGHS.Optimizer, env=nothing)
    # Create optimization model
    opt = isnothing(env) ? optimizer() : optimizer(env)
    model = direct_model(opt)
    set_silent(model)
    set_string_names_on_creation(model, false)
    # set_optimizer_attribute(model, "Threads", 4)

    # Define variables
    @variable(model, x[1:vertex_num])
    @variable(model, C)

    # Define constraints
    ϵ = 1  # Separation constant

    # Target MISs must have equal energy
    @constraint(model, target_set * x .== C)

    # Wrong MISs must have higher energy
    @constraint(model, wrong_set * x .<= C - ϵ)

    # All weights must be positive
    @constraint(model, [i=1:vertex_num], 1 <= x[i])

    # Solve the model
    optimize!(model)

    # Return results if feasible
    if is_solved_and_feasible(model)
        @info "Optimization successful!"
        return [value(x[i]) for i in 1:vertex_num]
    end

    # No feasible solution found
    return Float64[]
end

function _find_weight_new(mis_matrix::AbstractMatrix, prefix_buckets, optimizer=HiGHS.Optimizer, env=nothing)
    opt = isnothing(env) ? optimizer() : optimizer(env)
    model = direct_model(opt)
    set_silent(model)
    set_string_names_on_creation(model, false)
    # set_optimizer_attribute(model, "Threads", 4)

    y_dim, vertex_num = size(mis_matrix)

    @variable(model, x[1:vertex_num] >= 1, Int)
    @variable(model, y[1:y_dim], Bin)

    selected_indices = reduce(vcat, prefix_buckets)
    selected_set = Set(selected_indices)

    for idxs in prefix_buckets
        @constraint(model, sum(y[i] for i in idxs) == 1)
    end

    @constraint(model, [i=1:y_dim; i ∉ selected_set], y[i] == 0)

    row_sums = [sum(mis_matrix[i,:]) for i in 1:y_dim]
    M_big = 2 * maximum(row_sums)

    # Find the first selected index to use as a reference
    reference_indices = [first(idxs) for idxs in prefix_buckets]

    # Make all rows where y=1 equal to each other
    for i in selected_indices
        for ref_idx in reference_indices
            row_i = mis_matrix[i, :]
            row_ref = mis_matrix[ref_idx, :]
            # When both y_i and y_ref are 1, force their dot products to be equal
            # When either is 0, the constraint is relaxed
            @constraint(model, dot(row_i, x) - dot(row_ref, x) <= M_big * (2 - y[i] - y[ref_idx]))
            @constraint(model, dot(row_ref, x) - dot(row_i, x) <= M_big * (2 - y[i] - y[ref_idx]))
        end
    end

    # Make all rows where y=0 have lower values than rows where y=1
    for i in 1:y_dim
        for j in selected_indices
            row_i = mis_matrix[i, :]
            row_j = mis_matrix[j, :]
            # When y_i=0 and y_j=1, force row_i's dot product to be less than row_j's
            # Otherwise, the constraint is relaxed
            @constraint(model, dot(row_i, x) <= dot(row_j, x) - 1 + M_big * (1 - y[j] + y[i]))
        end
    end

    @objective(model, Min, sum(x))
    set_time_limit_sec(model, 30)
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

"""
    save_results_to_json(results::Vector{Gadget}, file_path::String)

Save search results to a JSON file.

# Arguments
- `results::Vector{Gadget}`: Vector of Gadget objects to save
- `file_path::String`: Path to save the JSON file

# Returns
- The file path where results were saved
"""
function save_results_to_json(results::Vector{<:AbstractGadget}, file_path::String)
    # Convert each result to a serializable dictionary
    json_results = map(results) do res
        base_dict = Dict(
            "rule_id" => res.rule_id,
            "ground_states" => res.ground_states,
            "pins" => res.pins,
            "graph" => Dict(
                "nodes" => [Dict("id" => i, "weight" => res.weights[i]) for i in 1:length(res.weights)],
                "edges" => [Dict("source" => src(e), "target" => dst(e)) for e in Graphs.edges(res.graph)],
                "graph_id" => res.graph_id
            )
        )

        # Add position data if it's a grid_gadget
        if isa(res, grid_gadget)
            base_dict["graph"]["positions"] = [Dict("id" => i, "position" => [pos[1], pos[2]]) for (i, pos) in enumerate(res.pos)]
        end

        return base_dict
    end

    # Write to JSON file with pretty formatting
    open(file_path, "w") do io
        write(io, JSON3.write(json_results; pretty=true))
    end

    return file_path
end
