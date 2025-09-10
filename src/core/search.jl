# Constants
const MAX_SUPPORTED_VERTICES = 32
const DEFAULT_EPSILON = 1.0
const DEFAULT_MAX_SAMPLES = 100

# Performance optimization: Cache for maximal independent sets
const MIS_CACHE = Dict{UInt64, Tuple{Vector{UInt32}, Int}}()

"""
    clear_cache!()

Clear the MIS cache to free memory. Useful for long-running processes.
"""
function clear_cache!()
    empty!(MIS_CACHE)
    GC.gc()  # Force garbage collection
    @info "MIS cache cleared, $(length(MIS_CACHE)) entries remaining"
end

"""
    get_cache_stats()

Get statistics about the MIS cache usage.
"""
function get_cache_stats()
    return (size=length(MIS_CACHE), memory_mb=Base.summarysize(MIS_CACHE) / 1024^2)
end

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
- `check_connectivity::Bool=true`: Whether to check graph connectivity after removing zero-weight vertices

# Returns
- `Tuple{Vector{Vector{Gadget}}, Vector{BitMatrix}}`: Found gadgets grouped by truth table and failed truth tables
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
    pin_candidates::Union{Nothing, Vector{Vector{Int}}}=nothing,
    check_connectivity::Bool=true
)
    # Initialize results as vector of vectors, one for each truth table
    results = Vector{Vector{Gadget}}(undef, length(truth_tables))
    for i in 1:length(truth_tables)
        results[i] = Gadget[]
    end
    failed_tt = BitMatrix[]
    total_found = 0

    # Performance monitoring
    start_time = time()
    initial_cache_size = length(MIS_CACHE)

    for (i, tt) in enumerate(truth_tables)
        tt_start = time()
        @info "searching for truth table: index $(i-1) [$total_found/$max_result_num found total]"
        
        # Calculate how many more results we need globally
        remaining_results = max_result_num - total_found
        if remaining_results <= 0
            @info "Max result number reached, stopping search"
            break
        end
        
        filter_fn = make_filter(tt, optimizer, env;
                                connected=connected,
                                objective=objective,
                                allow_defect=allow_defect,
                                max_samples=max_samples,
                                pin_candidates=pin_candidates,
                                check_connectivity=check_connectivity)
        # Pass remaining_results as max_results for early termination
        gadgets = find_matching_gadget(loader; filter=filter_fn, limit=limit, max_results=remaining_results)
        
        tt_time = time() - tt_start
        @info "Truth table $(i-1) processed in $(round(tt_time, digits=2))s, found $(length(gadgets)) gadgets"
        
        if !isempty(gadgets)
            # Store gadgets for this specific truth table
            results[i] = gadgets
            total_found += length(gadgets)
            
            # Save all results in flattened format for backward compatibility
            all_gadgets = vcat(results...)
            save_results_to_json(all_gadgets, save_path)
            
            # Break if we've reached the limit
            total_found >= max_result_num && break
        else
            push!(failed_tt, tt)
        end
        
        # Memory management: clear cache periodically for long searches
        if length(MIS_CACHE) > initial_cache_size + 1000
            @info "Cache growing large ($(length(MIS_CACHE)) entries), consider clearing if memory is an issue"
        end
    end

    total_time = time() - start_time
    cache_hits = length(MIS_CACHE) - initial_cache_size
    @info "Search completed in $(round(total_time, digits=2))s. Cache gained $cache_hits entries. Final cache stats: $(get_cache_stats())"

    return results, failed_tt
end


"""
    find_matching_gadget(loader; filter=nothing, limit=nothing, keys_range=nothing, max_results=nothing)

Find gadgets matching a given filter function from the graph loader.

# Arguments
- `loader::GraphLoader`: The graph loader containing candidate graphs

# Keywords
- `filter=nothing`: Filter function to apply to each graph
- `limit::Union{Int,Nothing}=nothing`: Maximum number of graphs to check
- `keys_range::Union{Nothing, Vector{Int}}=nothing`: Specific range of keys to search
- `max_results::Union{Int,Nothing}=nothing`: Maximum number of results to return (early termination)

# Returns
- `Vector{Gadget}`: Vector of matching gadgets
"""
function find_matching_gadget(loader::GraphLoader; filter=nothing, limit::Union{Int,Nothing}=nothing, keys_range::Union{Nothing, Vector{Int}}=nothing, max_results::Union{Int,Nothing}=nothing)
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
                # Early termination when max_results is reached
                if max_results !== nothing && length(results) >= max_results
                    @info "Early termination: found $(length(results)) results (max_results=$max_results)"
                    break
                end
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
- `check_connectivity::Bool=true`: Whether to check graph connectivity after removing zero-weight vertices

# Returns
- `Function`: Filter function that takes (graph, pos, pin_set) and returns (weights, truth_table, pin)
"""
function make_filter(truth_table::BitMatrix, optimizer, env; connected::Bool=false, objective=nothing, allow_defect::Bool=false, max_samples::Int=1000, pin_candidates::Union{Nothing, Vector{Vector{Int}}}=nothing, check_connectivity::Bool=true)
  
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
            weights = solve_weight_enumerate(maximal_independent_sets, target_mis_indices_all, vertex_num, candidate, optimizer, env, objective, allow_defect, max_samples, graph, check_connectivity)
            if !isempty(weights)
                return weights, truth_table, candidate
            end
        end
        return nothing, truth_table, nothing
    end
end


"""
    _graph_hash(g::SimpleGraph{Int}) -> UInt64

Create a hash for graph structure for caching purposes.
"""
function _graph_hash(g::SimpleGraph{Int})::UInt64
    # Create hash based on edges - more robust than adjacency matrix
    edges = sort!([(min(src(e), dst(e)), max(src(e), dst(e))) for e in Graphs.edges(g)])
    return hash((Graphs.nv(g), edges))
end

"""
    find_maximal_independent_sets(g)

Find all maximal independent sets of a graph using bit masks with caching.

# Arguments
- `g::SimpleGraph{Int}`: The input graph

# Returns
- `Tuple{Vector{UInt32}, Int}`: Bit masks representing maximal independent sets and their count

# Note
This function supports graphs with at most $(MAX_SUPPORTED_VERTICES) vertices due to UInt32 limitations.
Uses caching to avoid recomputing MIS for identical graphs.
"""
function find_maximal_independent_sets(g::SimpleGraph{Int})
    vertex_count = Graphs.nv(g)
    
    # Check vertex count limitation (UInt32 supports 32 vertices)
    if vertex_count > MAX_SUPPORTED_VERTICES
        error("Graph has $(vertex_count) vertices, but maximum supported is $(MAX_SUPPORTED_VERTICES). Consider reducing graph size or using matrix-based algorithms.")
    end
    
    # Check cache first
    graph_hash = _graph_hash(g)
    if haskey(MIS_CACHE, graph_hash)
        return MIS_CACHE[graph_hash]
    end
    
    # Compute complement once and reuse
    complement_g = Graphs.complement(g)
    cliques = Graphs.maximal_cliques(complement_g)
    
    # Pre-allocate with exact size
    masks = Vector{UInt32}(undef, length(cliques))
    
    # Vectorized bit operations for better performance
    @inbounds for (i, clique) in enumerate(cliques)
        mask::UInt32 = 0
        for v in clique
            mask |= UInt32(1) << (v - 1)
        end
        masks[i] = mask
    end
    
    result = (masks, length(masks))
    # Cache the result
    MIS_CACHE[graph_hash] = result
    
    return result
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
- `max_samples::Int=1000`: Maximum samples for enumeration
- `graph::Union{SimpleGraph{Int}, Nothing}=nothing`: Original graph for connectivity checking
- `check_connectivity::Bool=true`: Whether to check graph connectivity after removing zero-weight vertices

# Returns
- `Vector{Float64}`: Vertex weights if solution found, empty vector otherwise
"""
function solve_weight_enumerate(
    mis_result::Vector{UInt32},
    target_mis_indices_all::Vector{Vector{Int}},
    vertex_num::Int,
    pin_set::Vector{Int},
    optimizer,
    env=nothing,
    objective=nothing,
    allow_defect::Bool=false,
    max_samples::Int=1000,
    graph::Union{SimpleGraph{Int}, Nothing}=nothing,
    check_connectivity::Bool=true
)
    if optimizer === nothing
        error("Optimizer must be provided.")
    end

    any(isempty, target_mis_indices_all) && return Float64[]
    
    # Pre-compute constants to avoid repeated calculations
    all_mis_indices = collect(1:length(mis_result))
    
    # Early termination: estimate combination count and apply limit
    total_combinations = prod(length, target_mis_indices_all)
    if total_combinations > max_samples
        @warn "Too many combinations ($total_combinations), sampling first $max_samples"
        # Use iterative sampling instead of full enumeration
        return _solve_weight_sampled(mis_result, target_mis_indices_all, vertex_num, pin_set, 
                                   optimizer, env, objective, allow_defect, max_samples, all_mis_indices, graph, check_connectivity)
    end
    
    # Pre-allocate reusable containers
    target_indices_set = Vector{Int}(undef, length(target_mis_indices_all))
    target_set = Vector{UInt32}()
    wrong_set = Vector{UInt32}()
    
    # Enumerate all combinations to find weights
    sample_count = 0
    for target_indices in Iterators.product(target_mis_indices_all...)
        sample_count += 1
        if sample_count > max_samples
            break
        end
        
        # Reuse pre-allocated containers
        copyto!(target_indices_set, target_indices)
        wrong_indices = setdiff(all_mis_indices, target_indices_set)
        
        # Efficient array access
        resize!(target_set, length(target_indices_set))
        resize!(wrong_set, length(wrong_indices))
        
        @inbounds for (i, idx) in enumerate(target_indices_set)
            target_set[i] = mis_result[idx]
        end
        
        @inbounds for (i, idx) in enumerate(wrong_indices)
            wrong_set[i] = mis_result[idx]
        end

        weights = _find_weight(vertex_num, pin_set, target_set, wrong_set, optimizer, env, objective, allow_defect, graph, check_connectivity)
        if !isempty(weights)
            return weights
        end
    end

    return Float64[]
end

"""
    _solve_weight_sampled(...)

Optimized version for large combination spaces using random sampling.
"""
function _solve_weight_sampled(
    mis_result::Vector{UInt32},
    target_mis_indices_all::Vector{Vector{Int}},
    vertex_num::Int,
    pin_set::Vector{Int},
    optimizer,
    env,
    objective,
    allow_defect::Bool,
    max_samples::Int,
    all_mis_indices::Vector{Int},
    graph::Union{SimpleGraph{Int}, Nothing}=nothing,
    check_connectivity::Bool=true
)
    # Pre-allocate containers
    target_indices_set = Vector{Int}(undef, length(target_mis_indices_all))
    target_set = Vector{UInt32}()
    wrong_set = Vector{UInt32}()
    
    for _ in 1:max_samples
        # Random sampling from each dimension
        for (i, indices) in enumerate(target_mis_indices_all)
            target_indices_set[i] = rand(indices)
        end
        
        wrong_indices = setdiff(all_mis_indices, target_indices_set)
        
        # Efficient array construction
        resize!(target_set, length(target_indices_set))
        resize!(wrong_set, length(wrong_indices))
        
        @inbounds for (i, idx) in enumerate(target_indices_set)
            target_set[i] = mis_result[idx]
        end
        
        @inbounds for (i, idx) in enumerate(wrong_indices)
            wrong_set[i] = mis_result[idx]
        end

        weights = _find_weight(vertex_num, pin_set, target_set, wrong_set, optimizer, env, objective, allow_defect, graph, check_connectivity)
        if !isempty(weights)
            return weights
        end
    end
    
    return Float64[]
end

"""
    check_connectivity_after_removal(graph, vertices_to_remove)

Check if a graph remains connected after removing specified vertices.

# Arguments
- `graph::SimpleGraph{Int}`: The original graph
- `vertices_to_remove::Vector{Int}`: Vertices to remove (0-based indexing for weight array)

# Returns
- `Bool`: true if the graph remains connected after removal, false otherwise
"""
function check_connectivity_after_removal(graph::SimpleGraph{Int}, vertices_to_remove::Vector{Int})
    if isempty(vertices_to_remove)
        return Graphs.is_connected(graph)
    end
    
    # Convert to 1-based indexing for graph vertices
    vertices_to_remove_1based = vertices_to_remove .+ 1
    
    # Get remaining vertices
    all_vertices = collect(1:Graphs.nv(graph))
    remaining_vertices = setdiff(all_vertices, vertices_to_remove_1based)
    
    # If no vertices remain, it's not connected
    if isempty(remaining_vertices)
        return false
    end
    
    # Create subgraph with remaining vertices
    subgraph, _ = Graphs.induced_subgraph(graph, remaining_vertices)
    
    # Check if subgraph is connected
    return Graphs.is_connected(subgraph)
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
function _find_weight(vertex_num::Int, pin_set::Vector{Int}, target_masks::Vector{UInt32}, wrong_masks::Vector{UInt32}, optimizer, env, objective, allow_defect::Bool, graph::Union{SimpleGraph{Int}, Nothing}=nothing, check_connectivity::Bool=true)
    opt = isnothing(env) ? optimizer() : optimizer(env)
    model = direct_model(opt)
    set_silent(model)
    set_string_names_on_creation(model, false)

    if allow_defect
        x = @variable(model, [i=1:vertex_num], lower_bound = i in pin_set ? 1 : 0, Int)
    else
        @variable(model, x[1:vertex_num] >= 1, Int)
    end
    # @variable(model, x[1:vertex_num] >= 0)

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
        weights = [value(x[v]) for v in 1:vertex_num]
        
        # Check connectivity after removing zero-weight vertices (if enabled)
        if check_connectivity && graph !== nothing
            # Find vertices with zero (or near-zero) weights
            zero_weight_vertices = findall(w -> abs(w) < 1e-6, weights)
            
            if !isempty(zero_weight_vertices)
                # Convert to 0-based indexing for connectivity check
                zero_vertices_0based = zero_weight_vertices .- 1
                
                if !check_connectivity_after_removal(graph, zero_vertices_0based)
                    @info "Solution found but graph becomes disconnected after removing zero-weight vertices, rejecting"
                    return Float64[]
                end
            end
        end
        
        @info "found a valid solution with connectivity preserved"
        @show weights
        return weights
    end
    return Float64[]
end