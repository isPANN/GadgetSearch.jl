# Constants
const MAX_SUPPORTED_VERTICES = 32
const DEFAULT_EPSILON = 1.0
const DEFAULT_MAX_SAMPLES = 100

"""
    search_by_truth_tables(loader, truth_tables; kwargs...)

Search for multiple truth tables by reusing the graph search functionality.

# Arguments
- `loader::GraphLoader`: The graph loader containing candidate graphs
- `truth_tables::Vector{BitMatrix}`: Truth tables to search for

# Keywords
- `bit_num`: Number of bits for the gadget
- `optimizer`: Optimization solver to use for weight finding
- `env=nothing`: Environment for the optimizer
- `connected::Bool=false`: Whether to require connected graphs only
- `objective=nothing`: Objective function for optimization
- `allow_defect::Bool=false`: Whether to allow defective gadgets
- `limit=nothing`: Maximum number of graphs to search
- `max_samples::Int=100`: Maximum samples for weight enumeration
- `save_path::String="results.json"`: Path to save intermediate results
- `pin_candidates::Union{Nothing, Vector{Vector{Int}}}=nothing`: Candidate pin combinations

# Returns
- `Tuple{Vector{Gadget}, Vector{BitMatrix}}`: Found gadgets and failed truth tables
"""
function search_by_truth_tables(
    loader::GraphLoader,
    truth_tables::Vector{BitMatrix};
    optimizer,
    env=nothing,
    connected::Bool=false,
    objective=nothing,
    allow_defect::Bool=false,
    limit=nothing,
    max_samples::Int=DEFAULT_MAX_SAMPLES,
    max_result_num::Int=1,
    save_path::String="results.json",
    pin_candidates::Union{Nothing, Vector{Vector{Int}}}=nothing
)
    results = Gadget[]
    failed_tt = BitMatrix[]

    for (i, tt) in enumerate(truth_tables)
        @info "searching for truth table: index $(i-1)"
        filter_fn = make_filter(tt, optimizer, env;
                                connected=connected,
                                objective=objective,
                                allow_defect=allow_defect,
                                max_samples=max_samples,
                                pin_candidates=pin_candidates)
        gadgets = find_matching_gadget(loader; filter=filter_fn, limit=limit)
        if !isempty(gadgets)
            # Add gadgets but respect the max_result_num limit
            for gadget in gadgets
                push!(results, gadget)
                if length(results) >= max_result_num
                    break
                end
            end
            save_results_to_json(results, save_path)
            # Break if we've reached the limit
            length(results) >= max_result_num && break
        else
            push!(failed_tt, tt)
        end
    end

    return results, failed_tt
end


"""
    find_matching_gadget(loader; filter=nothing, limit=nothing, keys_range=nothing)

Find gadgets matching a given filter function from the graph loader.

# Arguments
- `loader::GraphLoader`: The graph loader containing candidate graphs

# Keywords
- `filter=nothing`: Filter function to apply to each graph
- `limit::Union{Int,Nothing}=nothing`: Maximum number of graphs to check
- `keys_range::Union{Nothing, Vector{Int}}=nothing`: Specific range of keys to search

# Returns
- `Vector{Gadget}`: Vector of matching gadgets
"""
function find_matching_gadget(loader::GraphLoader; filter=nothing, limit::Union{Int,Nothing}=nothing, keys_range::Union{Nothing, Vector{Int}}=nothing)
    keys_raw = keys_range === nothing ? collect(keys(loader)) : keys_range
    keys_to_search = isa(keys_raw[1], Int) ? keys_raw : parse.(Int, keys_raw)
    total = limit === nothing ? length(keys_to_search) : min(length(keys_to_search), limit)

    results = Gadget[]

    @showprogress for key in Iterators.take(keys_to_search, total)
        g = loader[key]
        if filter !== nothing
            weights, truth_table, pin = filter(g, loader.layout[key], loader.pinset)
            if !isnothing(weights)
                push!(results, Gadget(truth_table, g, pin, weights, loader.layout[key])) 
            end
        end
    end
    return results
end

"""
    make_filter(truth_table, optimizer, env; kwargs...)

Create a filter closure for graph search that checks if a graph can implement the given truth table.

# Arguments
- `truth_table::BitMatrix`: The target truth table to implement
- `optimizer`: Optimization solver for weight finding
- `env`: Environment for the optimizer

# Keywords
- `connected::Bool=false`: Whether to require connected graphs only
- `objective=nothing`: Objective function for optimization
- `allow_defect::Bool=false`: Whether to allow defective gadgets
- `max_samples::Int=1000`: Maximum samples for weight enumeration
- `pin_candidates::Union{Nothing, Vector{Vector{Int}}}=nothing`: Candidate pin combinations

# Returns
- `Function`: Filter function that takes (graph, pos, pin_set) and returns (weights, truth_table, pin)
"""
function make_filter(truth_table::BitMatrix, optimizer, env; connected::Bool=false, objective=nothing, allow_defect::Bool=false, max_samples::Int=1000, pin_candidates::Union{Nothing, Vector{Vector{Int}}}=nothing)
  
    return function(graph::SimpleGraph{Int}, pos::Vector{Tuple{Float64, Float64}}, pin_set::Union{Nothing, Vector{Int}}=nothing)
        if !connected
            Graphs.is_connected(graph) || return nothing, truth_table, nothing
        end

        vertex_num = Graphs.nv(graph)

        maximal_independent_sets, set_count = find_maximal_independent_sets(graph)

        if pin_set === nothing
            pin_set = collect(1:vertex_num)
        end
        
        if pin_candidates !== nothing
            all_candidates = Vector{Vector{Int}}()
            for candidate in pin_candidates
                if all(pin in pin_set for pin in candidate)
                    push!(all_candidates, candidate)
                else
                    @warn "Candidate $candidate contains pins not in pin_set $pin_set, skipping"
                end
            end
        else
            @warn "pin_candidates is not provided, using all possible pin combinations"
            all_candidates = collect(Combinatorics.combinations(pin_set, size(truth_table, 2)))
        end
        
        for candidate in all_candidates
            target_mis_indices_all = match_rows_by_pinset(maximal_independent_sets, truth_table, candidate)
            weights = solve_weight_enumerate(maximal_independent_sets, target_mis_indices_all, vertex_num, optimizer, env, objective, allow_defect, max_samples)
            if !isempty(weights)
                return weights, truth_table, candidate
            end
        end
        return nothing, truth_table, nothing
    end
end


"""
    find_maximal_independent_sets(g)

Find all maximal independent sets of a graph using bit masks.

# Arguments
- `g::SimpleGraph{Int}`: The input graph

# Returns
- `Tuple{Vector{UInt32}, Int}`: Bit masks representing maximal independent sets and their count

# Note
This function supports graphs with at most $(MAX_SUPPORTED_VERTICES) vertices due to UInt32 limitations.
"""
function find_maximal_independent_sets(g::SimpleGraph{Int})
    cliques = Graphs.maximal_cliques(Graphs.complement(g))
    vertex_count = Graphs.nv(g)
    
    # Check vertex count limitation (UInt32 supports 32 vertices)
    if vertex_count > MAX_SUPPORTED_VERTICES
        error("Graph has $(vertex_count) vertices, but maximum supported is $(MAX_SUPPORTED_VERTICES). Consider reducing graph size or using matrix-based algorithms.")
    end
    
    # Use UInt32 for efficient bit operations
    masks = UInt32[]
    sizehint!(masks, length(cliques))  # Pre-allocate memory
    
    for clique in cliques
        mask::UInt32 = 0
        for v in clique
            mask |= UInt32(1) << (v - 1)
        end
        push!(masks, mask)
    end
    
    return masks, length(masks)
end

"""
    match_rows_by_pinset(masks, truth_table, pin_set)

Match truth table rows to maximal independent sets based on pin configurations.

# Arguments
- `masks::Vector{UInt32}`: Bit masks representing maximal independent sets
- `truth_table::BitMatrix`: The target truth table
- `pin_set::Vector{Int}`: Pin positions to consider

# Returns
- `Vector{Vector{Int}}`: For each truth table row, indices of matching maximal independent sets
"""
function match_rows_by_pinset(masks::Vector{UInt32}, truth_table::BitMatrix, pin_set::Vector{Int})::Vector{Vector{Int}}
    num_rows = size(truth_table, 1)
    result = Vector{Vector{Int}}(undef, num_rows)

    for i in 1:num_rows
        # Build query mask
        query_mask::UInt32 = 0
        for (bit_pos, pin) in enumerate(pin_set)
            query_mask |= UInt32(truth_table[i, bit_pos]) << (bit_pos - 1)
        end
        
        # Find matching MIS
        matches = Int[]
        for (j, m) in enumerate(masks)
            # Extract values at pin positions from MIS
            extracted::UInt32 = 0
            for (bit_pos, pin) in enumerate(pin_set)
                extracted |= ((m >> (pin - 1)) & 0x1) << (bit_pos - 1)
            end
            
            if extracted == query_mask
                push!(matches, j)
            end
        end
        result[i] = matches
    end

    return result
end

"""
    solve_weight_enumerate(mis_result, target_mis_indices_all, vertex_num, optimizer; kwargs...)

Solve for vertex weights by enumerating all combinations of target MIS indices.

# Arguments
- `mis_result::Vector{UInt32}`: Maximal independent sets as bit masks
- `target_mis_indices_all::Vector{Vector{Int}}`: Target MIS indices for each truth table row
- `vertex_num::Int`: Number of vertices in the graph
- `optimizer`: Optimization solver to use

# Keywords
- `env=nothing`: Environment for the optimizer
- `objective=nothing`: Objective function for optimization
- `allow_defect::Bool=false`: Whether to allow defective gadgets
- `max_samples::Int=1000`: Maximum samples for enumeration (currently unused)

# Returns
- `Vector{Float64}`: Vertex weights if solution found, empty vector otherwise
"""
function solve_weight_enumerate(
    mis_result::Vector{UInt32},
    target_mis_indices_all::Vector{Vector{Int}},
    vertex_num::Int,
    optimizer,
    env=nothing,
    objective=nothing,
    allow_defect::Bool=false,
    max_samples::Int=1000
)
    if optimizer === nothing
        error("Optimizer must be provided.")
    end

    any(isempty, target_mis_indices_all) && return Float64[]
    all_mis_indices = eachindex(mis_result)
    
    # Enumerate all combinations to find weights
    for target_indices in Iterators.product(target_mis_indices_all...)
        target_indices_set = collect(target_indices)
        wrong_indices = setdiff(all_mis_indices, target_indices_set)
        target_set = mis_result[target_indices_set]
        wrong_set = mis_result[wrong_indices]

        weights = _find_weight(vertex_num, target_set, wrong_set, optimizer, env, objective, allow_defect)
        if !isempty(weights)
            return weights
        end
    end

    return Float64[]
end

"""
    _find_weight(vertex_num, target_masks, wrong_masks, optimizer, env, objective, allow_defect)

Find vertex weights using optimization to distinguish target and wrong maximal independent sets.

# Arguments
- `vertex_num::Int`: Number of vertices
- `target_masks::Vector{UInt32}`: Bit masks for target MIS (should have equal energy)
- `wrong_masks::Vector{UInt32}`: Bit masks for wrong MIS (should have higher energy)
- `optimizer`: Optimization solver
- `env`: Environment for the optimizer
- `objective`: Objective function
- `allow_defect::Bool`: Whether to allow defective gadgets

# Returns
- `Vector{Float64}`: Vertex weights if solution found, empty vector otherwise
"""
function _find_weight(vertex_num::Int, target_masks::Vector{UInt32}, wrong_masks::Vector{UInt32}, optimizer, env, objective, allow_defect::Bool)
    opt = isnothing(env) ? optimizer() : optimizer(env)
    model = direct_model(opt)
    set_silent(model)
    set_string_names_on_creation(model, false)

    # if allow_defect
    #     x = @variable(model, [i=1:vertex_num], lower_bound = i <= 4 ? 1 : 0)
    # else
    #     @variable(model, x[1:vertex_num] >= 1)
    # end
    @variable(model, x[1:vertex_num] >= 1)

    @variable(model, C)

    # Target MIS must have equal energy
    for m in target_masks
        @constraint(model, sum(((m >> (v - 1)) & 0x1) * x[v] for v in 1:vertex_num) == C)
    end

    # Wrong MIS must have higher energy
    for m in wrong_masks
        @constraint(model, sum(((m >> (v - 1)) & 0x1) * x[v] for v in 1:vertex_num) <= C - DEFAULT_EPSILON)
    end

    if objective !== nothing
        @objective(model, Min, objective(x))
    end

    optimize!(model)
    if is_solved_and_feasible(model)
        @info "found a solution"
        return [value(x[v]) for v in 1:vertex_num]
    end
    return Float64[]
end