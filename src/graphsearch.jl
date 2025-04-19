# Structure to encapsulate search parameters for better organization and type safety
"""
    SearchParameters(;
        pin_set::Vector{Int} = Int[],
        max_file_size_mb::Int = 30,
        split_size::Int = 700_000,
        start_idx::Int = 0,
        end_idx::Int = 0,
        greedy::Bool = false,
        threshold::Int = 0,
        max_samples::Int = 0
    )

Contains parameters for graph search operations

# Fields
- `pin_set::Vector{Int}`: Set of pins to search for
- `max_file_size_mb::Int`: Maximum file size in MB to process before splitting
- `split_size::Int`: Maximum number of rows in each split file
- `start_idx::Int`: Starting index for graph search
- `end_idx::Int`: Ending index for graph search
- `greedy::Bool`: Whether to use greedy approach
- `threshold::Int`: Maximum number of MIS combinations to consider
- `max_samples::Int`: Maximum number of samples to generate
"""
struct SearchParameters
    pin_set::Vector{Int}
    max_file_size_mb::Int
    split_size::Int
    start_idx::Int
    end_idx::Int
    greedy::Bool
    threshold::Int
    max_samples::Int

    # Constructor with default values and validation
    function SearchParameters(;
        pin_set::Vector{Int} = Int[],
        max_file_size_mb::Int = 30,
        split_size::Int = 700_000,
        start_idx::Int = 0,
        end_idx::Int = 0,
        greedy::Bool = false,
        threshold::Int = 0,
        max_samples::Int = 0
    )
        # Parameter validation
        max_file_size_mb >= 0 || throw(ArgumentError("max_file_size_mb must be positive"))
        split_size > 0 || throw(ArgumentError("split_size must be positive"))
        start_idx >= 0 || throw(ArgumentError("start_idx must be non-negative"))
        end_idx >= 0 || throw(ArgumentError("end_idx must be non-negative"))
        threshold >= 0 || throw(ArgumentError("threshold must be non-negative"))
        max_samples >= 0 || throw(ArgumentError("max_samples must be non-negative"))

        new(pin_set, max_file_size_mb, split_size, start_idx, end_idx, greedy, threshold, max_samples)
    end
end

"""
    search_rules(graph_path, bit_num, gate_list, save_path;
                    pin_set::Vector{Int} = Int[],
                    # File processing parameters
                    max_file_size_mb::Int=30, split_size::Int=700_000,
                    # Search strategy parameters
                    start_idx::Int=0, end_idx::Int=0,
                    greedy::Bool=false, threshold::Int=0, max_samples::Int=0,
                    optimizer=HiGHS.Optimizer)

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
                    optimizer=HiGHS.Optimizer)
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

    # Track rules that couldn't be found
    non_results = Int[]
    save_res = Vector{Gadget}()

    # Process each rule in the gate list
    for rule_id in gate_list
        # Search for a single rule with the given parameters
        res = search_single_rule(graph_path, bit_num;
                               rule_id=rule_id,
                               optimizer=optimizer,
                               pin_set=params.pin_set,
                               max_file_size_mb=params.max_file_size_mb,
                               split_size=params.split_size,
                               start_idx=params.start_idx,
                               end_idx=params.end_idx,
                               greedy=params.greedy,
                               threshold=params.threshold,
                               max_samples=params.max_samples)

        if !isnothing(res)
            push!(save_res, res)
            save_results_to_json(save_res, save_path)
        else
            push!(non_results, rule_id)
        end
    end

    return non_results
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
                            optimizer=HiGHS.Optimizer)::Union{Gadget, Nothing}

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
                            optimizer=HiGHS.Optimizer
)::Union{Gadget, Nothing}
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
    isempty(ground_states) && return nothing

    # Check if file needs to be split
    if (filesize(graph_path) / (1024 * 1024)) > params.max_file_size_mb
        # Split the file into smaller parts and process each part separately
        split_file_paths = _split_large_file(graph_path, params.split_size)

        # Process each split file
        for (i, file) in enumerate(split_file_paths)
            graph_dict = read_graph_dict(file)

            # Execute search on this part of the graph database
            result = _execute_graph_search(
                graph_dict, bit_num, ground_states, rule_id,
                params.pin_set, i, params.split_size, params.start_idx,
                params.end_idx, params.greedy, params.threshold,
                params.max_samples, optimizer
            )

            # Return early if a solution is found
            isnothing(result) || return result
        end
    else
        # Process the entire file at once
        graph_dict = read_graph_dict(graph_path)
        result = _execute_graph_search(
            graph_dict, bit_num, ground_states, rule_id,
            params.pin_set, 0, params.split_size, params.start_idx,
            params.end_idx, params.greedy, params.threshold,
            params.max_samples, optimizer
        )
        return result
    end

    # No solution found
    return nothing
end

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
                            optimizer=HiGHS.Optimizer)
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
        graph_dict, bit_num, ground_states, rule_id,
        params.pin_set, split_idx, params.split_size, params.start_idx,
        params.end_idx, params.greedy, params.threshold,
        params.max_samples, optimizer
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
                            optimizer=HiGHS.Optimizer)
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

    # Search on the single graph
    pins, weight = _search_on_single_graph(graph, bit_num, ground_states, params, optimizer)

    # Return result if found
    if !isnothing(pins)
        return convert_to_result(0, graph, pins, weight, ground_states, rule_id)
    else
        @info "No valid solution found in this graph."
        return nothing
    end
end


function _execute_graph_search(graph_dict::Dict{String, SimpleGraph{Int}},
                                bit_num::Union{Int, Vector{Int}},
                                ground_state::Vector{Int},
                                rule_id::Int=-1,
                                pin_set::Vector{Int}=Int[],
                                split_idx::Int=0, split_size::Int=700_000,
                                start_idx::Int=0, end_idx::Int=0,
                                greedy::Bool=false, threshold::Int=0, max_samples::Int=0,
                                optimizer=HiGHS.Optimizer)
    # Calculate base offset for graph IDs
    base_offset = split_idx == 0 ? 0 : (split_idx - 1) * split_size

    # Filter and enumerate the graphs to process
    filtered_graphs = Iterators.enumerate(graph_dict) |>
        iter -> Iterators.filter(((i, _),) ->
            (start_idx <= 0 || base_offset + i >= start_idx) &&
            (end_idx <= 0 || base_offset + i <= end_idx), iter)

    # Process each graph
    for (i, (gname, graph)) in filtered_graphs
        # Print progress message every 10,000 graphs
        (base_offset + i) % 10000 == 0 && println("Searched $(base_offset + i) graphs...")

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

        # Perform the search on the current graph
        pins, weight = _search_on_single_graph(graph, bit_num, ground_state, params, optimizer)

        # If a valid solution is found, return the result
        if !isnothing(pins)
            # Calculate graph ID based on split information
            graph_id = _extract_numbers(gname) + (split_idx == 0 ? 0 : base_offset)
            return convert_to_result(graph_id, graph, pins, weight, ground_state, rule_id)
        end
    end

    # No solution found
    return nothing
end

"""_search_on_single_graph(graph::SimpleGraph{Int}, bit_num::Int, ground_states::Vector{Int}, params::SearchParameters, optimizer=HiGHS.Optimizer)

Search for a solution on a single graph for generic constraints.

# Arguments
- `graph::SimpleGraph{Int}`: The graph to search
- `bit_num::Int`: Number of bits
- `ground_states::Vector{Int}`: Ground states to satisfy
- `params::SearchParameters`: Search parameters including pin_set, greedy, threshold, and max_samples
- `optimizer`: Optimization solver to use

# Returns
- Tuple of (pins, weights) if a solution is found, otherwise (nothing, nothing)
"""
function _search_on_single_graph(graph::SimpleGraph{Int}, bit_num::Int, ground_states::Vector{Int}, params::SearchParameters, optimizer=HiGHS.Optimizer)
    # This is a version for generic constraints.
    # Ground states form a vector by converting binary numbers (e.g. 101) to their corresponding decimal values (e.g. 5).
    @assert length(ground_states) > 0 && maximum(ground_states) < 2^bit_num

    # Graph must be connected.
    Graphs.is_connected(graph) || return nothing, nothing
    vertex_num = Graphs.nv(graph)

    # Find all maximal independent sets of this graph.
    mis_result, _ = find_maximal_independent_sets(graph)

    # Generate all pin vectors, i.e. the permutations of `bit_num` vertices in `pin_set`.
    local pin_set = params.pin_set
    if isempty(pin_set)
        pin_set = collect(1:vertex_num)
    end
    @assert length(pin_set) >= bit_num
    all_candidates = _generate_constraint_bit(pin_set, bit_num)

    # Iterate over all possible pin vectors and search the vertex weight.
    for candidate in all_candidates
        # Check if the `ground_states` are contained in the MISs under the choice of `candidate`.
        target_mis_indices_all = _check_grstates_in_candidate(mis_result, candidate, ground_states, params.greedy)
        isempty(target_mis_indices_all) && continue

        # Sample combinations if threshold and max_samples are set
        local selected_combinations
        if params.threshold > 0 && params.max_samples > 0
            selected_combinations = _sample_possible_mis(target_mis_indices_all, params.threshold, params.max_samples)
        else
            selected_combinations = collect(Iterators.product(target_mis_indices_all...))
        end

        # Try each combination
        for combination in selected_combinations
            # Split MIS set into target and wrong sets
            target_mis_set, wrong_mis_set = _split_mis_set(mis_result, vcat(combination...))

            # Find weights
            weight = _find_weight(vertex_num, target_mis_set, wrong_mis_set, optimizer)
            isempty(weight) && continue

            # Return solution if found
            return candidate, weight
        end
    end
    return nothing, nothing
end


"""_search_on_single_graph(graph::SimpleGraph{Int}, bit_num::Vector{Int}, ground_states::Vector{Int}, params::SearchParameters, optimizer=HiGHS.Optimizer)

Search for a solution on a single graph for logic gates.

# Arguments
- `graph::SimpleGraph{Int}`: The graph to search
- `bit_num::Vector{Int}`: Vector of [input_bits, output_bits]
- `ground_states::Vector{Int}`: Ground states to satisfy
- `params::SearchParameters`: Search parameters including pin_set, greedy, threshold, and max_samples
- `optimizer`: Optimization solver to use

# Returns
- Tuple of (pins, weights) if a solution is found, otherwise (nothing, nothing)
"""
function _search_on_single_graph(graph::SimpleGraph{Int}, bit_num::Vector{Int}, ground_states::Vector{Int}, params::SearchParameters, optimizer=HiGHS.Optimizer)
    # This is a version for a logic gate.
    @assert length(bit_num) == 2
    @assert length(ground_states) > 0 && maximum(ground_states) < 2^(sum(bit_num))
    input_num = bit_num[1]
    output_num = bit_num[2]

    # Graph must be connected
    Graphs.is_connected(graph) || return nothing, nothing
    vertex_num = Graphs.nv(graph)

    # Find all maximal independent sets of this graph
    mis_result, mis_num = find_maximal_independent_sets(graph)

    # Check if there are enough MISs for the input bits
    mis_num < 2^(input_num) && return nothing, nothing

    # Generate input candidates
    local pin_set = params.pin_set
    if isempty(pin_set)
        # If the pin set is not provided, all vertices are considered as potential pins
        pin_set = collect(1:vertex_num)
        # Generate valid input pin combinations
        input_candidates = _generate_pin_set(graph, mis_result, input_num)
    else
        # If the pin set is provided, all selections in `pin_set` are considered valid by default
        input_candidates = _generate_gate_input(pin_set, input_num)
    end

    # Check if there are any valid input candidates
    isempty(input_candidates) && return nothing, nothing

    # Try each input candidate
    for candidate in input_candidates
        # Find remaining elements for output pins
        remain_elements = setdiff(pin_set, candidate)

        # Try each possible output pin combination
        for output_bits in permutations(remain_elements, output_num)
            # Combine input and output pins
            candidate_full = vcat(candidate, output_bits)

            # Check if ground states are contained in MISs
            target_mis_indices_all = _check_grstates_in_candidate(mis_result, candidate_full, ground_states, params.greedy)
            isempty(target_mis_indices_all) && continue

            # Sample combinations if threshold and max_samples are set
            local selected_combinations
            if params.threshold > 0 && params.max_samples > 0
                selected_combinations = _sample_possible_mis(target_mis_indices_all, params.threshold, params.max_samples)
            else
                selected_combinations = collect(Iterators.product(target_mis_indices_all...))
            end

            # Try each combination
            for combination in selected_combinations
                # Split MIS set into target and wrong sets
                target_mis_set, wrong_mis_set = _split_mis_set(mis_result, vcat(combination...))

                # Find weights
                weight = _find_weight(vertex_num, target_mis_set, wrong_mis_set, optimizer)
                isempty(weight) && continue

                # Return solution if found
                return candidate_full, weight
            end
        end
    end
    return nothing, nothing
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
function _search_on_single_graph(graph::SimpleGraph{Int}, bit_num::Vector{Int}, rule_id::Int; pin_set::Vector{Int}=Int[], greedy::Bool=true, threshold::Int=0, max_samples::Int=0, optimizer=HiGHS.Optimizer)
    @assert length(bit_num) == 2
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
    return _search_on_single_graph(graph, bit_num, generic_rule(rule_id, bit_num), params, optimizer)
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

"""_check_grstates_in_candidate(mis_result::AbstractMatrix{Int}, candidate::Vector{Int}, grstates::Vector{Int}, greedy::Bool=false)::Vector{Vector{Int}}

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
function _find_weight(vertex_num::Int, target_set::AbstractMatrix{Int}, wrong_set::AbstractMatrix{Int}, optimizer=HiGHS.Optimizer)
    # Create optimization model
    model = direct_model(optimizer())
    set_silent(model)
    set_string_names_on_creation(model, false)

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


"""convert_to_result(graph_id::Int, g::SimpleGraph{Int}, pins::Vector{Int}, weight::Vector{T}, ground_states::Vector{Int}, rule_id::Int=0) where T <: Real

Convert search results to a Gadget object.

# Arguments
- `graph_id::Int`: ID of the graph
- `g::SimpleGraph{Int}`: The graph
- `pins::Vector{Int}`: Pin vertices
- `weight::Vector{T}`: Vertex weights
- `ground_states::Vector{Int}`: Ground states
- `rule_id::Int=0`: ID of the logic gate

# Returns
- A Gadget object containing the search results
"""
function convert_to_result(graph_id::Int, g::SimpleGraph{Int}, pins::Vector{Int}, weight::Vector{T}, ground_states::Vector{Int}, rule_id::Int=0) where T <: Real
    # Format ground states for output
    ground_state_io = format_grstate_output(ground_states, length(pins))

    # Format edges and nodes for display
    edges_str = join(["$(src(e)) -- $(dst(e))" for e in Graphs.edges(g)], ", ")
    nodes_str = join(["$i -- $(weight[i])" for i in 1:length(weight)], ", ")

    # Log the result
    @info """ === Result ===
    Gate ID: $rule_id
    Ground States: $ground_state_io
    Graph ID: $graph_id
    Edges: $edges_str
    Nodes: $nodes_str
    """

    # Create and return the Gadget object
    return Gadget(
        rule_id,
        ground_state_io,
        graph_id,
        g,
        pins,
        weight
    )
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
function save_results_to_json(results::Vector{Gadget}, file_path::String)
    # Convert each result to a serializable dictionary
    json_results = map(results) do res
        Dict(
            "rule_id" => res.rule_id,
            "ground_states" => res.ground_states,
            "pins" => res.pins,
            "graph" => Dict(
                "nodes" => [Dict("id" => i, "weight" => res.weights[i]) for i in 1:length(res.weights)],
                "edges" => [Dict("source" => src(e), "target" => dst(e)) for e in Graphs.edges(res.graph)],
                "graph_id" => res.graph_id
            )
        )
    end

    # Write to JSON file with pretty formatting
    open(file_path, "w") do io
        write(io, JSON3.write(json_results; pretty=true))
    end

    return file_path
end

