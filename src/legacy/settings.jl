# Trait system for GadgetSearch

"""
Abstract trait for search strategies.
"""
abstract type AbstractSearchStrategy end

"""
Trait for generic constraint search.
"""
struct GenericConstraintSearch <: AbstractSearchStrategy end

"""
Trait for logic gate search.
"""
struct LogicGateSearch <: AbstractSearchStrategy end

"""
Determine the appropriate gadget type based on graph type trait.
"""
function gadget_type(::GeneralGraph)
    return Gadget
end

function gadget_type(::GridGraph)
    return grid_gadget
end

"""
Create an empty result vector based on graph type trait.
"""
function create_result_vector(::GeneralGraph)
    return Vector{Gadget}()
end

function create_result_vector(::GridGraph)
    return Vector{grid_gadget}()
end

"""
Convert search results to a gadget object based on graph type trait.
"""
function convert_to_gadget(::GeneralGraph, graph_id::Int, g::SimpleGraph{Int}, pins::Vector{Int},
                          weight::Vector{T}, ground_states::Vector{Int}, rule_id::Int) where T <: Real
    # Format ground states for output
    ground_state_io = format_grstate_output(ground_states, length(pins))

    return Gadget(
        rule_id,
        ground_state_io,
        graph_id,
        g,
        pins,
        weight
    )
end

function convert_to_gadget(graph_type::GridGraph, graph_id::Int, g::SimpleGraph{Int}, pins::Vector{Int},
                          weight::Vector{T}, ground_states::Vector{Int}, rule_id::Int) where T <: Real
    # Format ground states for output
    ground_state_io = format_grstate_output(ground_states, length(pins))

    # Check if we have a path to position data
    if !isnothing(graph_type.pos_data_path) && isfile(graph_type.pos_data_path)
        # Load position data
        pos_data = JSON.parsefile(graph_type.pos_data_path)
        # Check if we have position data for this graph
        if haskey(pos_data, "$graph_id")
            # Get position data for this graph
            pos_int = Int.(pos_data["$graph_id"])
            # Convert position indices to tuples
            pos = _index_to_tuple(pos_int, graph_type.grid_dims)

            return grid_gadget(
                rule_id,
                ground_state_io,
                graph_id,
                g,
                pins,
                weight,
                pos
            )
        end
    end

    # Fallback to regular Gadget if position data is not available
    @warn "Position data not available for graph $graph_id, returning regular Gadget"
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
Determine the search strategy based on bit_num.
"""
function determine_search_strategy(bit_num::Int)
    return GenericConstraintSearch()
end

function determine_search_strategy(bit_num::Vector{Int})
    if length(bit_num) == 2
        return LogicGateSearch()
    else
        return GenericConstraintSearch()
    end
end

"""
Search on a single graph based on search strategy trait.
"""
function search_on_single_graph(::GenericConstraintSearch, graph::SimpleGraph{Int}, bit_num::Int,
                               ground_states::Vector{Int}, params::SearchParameters, optimizer, env)
    # Implementation for generic constraint search
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

        weight = _find_weight_new(mis_result, target_mis_indices_all, optimizer, env)
        isempty(weight) && continue
        return candidate, weight
    end
    return nothing, nothing
end

function search_on_single_graph(::LogicGateSearch, graph::SimpleGraph{Int}, bit_num::Vector{Int},
                               ground_states::Vector{Int}, params::SearchParameters, optimizer, env)
    # Implementation for logic gate search
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
                weight = _find_weight(vertex_num, target_mis_set, wrong_mis_set, optimizer, env)
                isempty(weight) && continue

                # Return solution if found
                return candidate_full, weight
            end

            # weight = _find_weight_new(mis_result, target_mis_indices_all, optimizer, env)
            # isempty(weight) && continue
            # return candidate_full, weight
        end
    end
    return nothing, nothing
end
