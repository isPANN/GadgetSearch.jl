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

# ============================================================================
# Abstract Types for Energy Models
# ============================================================================

"""
    EnergyModel

Abstract type for energy models. Subtypes define how energy is computed from vertex states.
"""
abstract type EnergyModel end

"""
    RydbergModel <: EnergyModel

Energy model for Rydberg atom systems.
- State space: Maximal Independent Sets (MIS)
- Energy: E(σ) = Σᵢ hᵢσᵢ (vertex weights only)
"""
struct RydbergModel <: EnergyModel end

"""
    QUBOModel <: EnergyModel

Energy model for general QUBO (Quadratic Unconstrained Binary Optimization).
- State space: All 2^n binary configurations
- Energy: E(σ) = Σᵢ hᵢσᵢ + Σ₍ᵢ,ⱼ₎ Jᵢⱼσᵢσⱼ (vertex + edge weights)
"""
struct QUBOModel <: EnergyModel end

"""
    RydbergUnweightedModel <: EnergyModel

Energy model for unweighted Rydberg atom systems using reduced α-tensor verification.

- State space: Maximal Independent Sets (MIS)
- Verification: Compute reduced α-tensor, determine ground state boundary configurations
  by maximizing total(s) = α[s] + count_ones(s), and check they match the truth table.
- No optimizer needed: pure combinatorial verification via graph structure
- All vertex weights are uniform (= 1)

The reduced α-tensor α[s] gives the maximum IS size in the interior of the graph
for each boundary (pin) configuration s. Ground states are the configurations
that maximize the total IS size (interior + boundary).

Reference: Liu et al., "Computer-assisted gadget design and problem reduction
of unweighted maximum independent set", PRX Quantum 4, 010316 (2023)
"""
struct RydbergUnweightedModel <: EnergyModel end

# ============================================================================
# Abstract Types for Constraints
# ============================================================================

"""
    GadgetConstraint

Abstract type for gadget constraints. Subtypes define what ground states are required.
"""
abstract type GadgetConstraint end

"""
    TruthTableConstraint <: GadgetConstraint

Constraint defined by a truth table (for Rydberg/MIS-based gadgets).

# Fields
- `truth_table::BitMatrix`: Truth table where each row is a ground state configuration on pins
"""
struct TruthTableConstraint <: GadgetConstraint
    truth_table::BitMatrix
    
    function TruthTableConstraint(tt::BitMatrix)
        new(tt)
    end
end

# Convenience constructor from BitMatrix
TruthTableConstraint(tt::Matrix{Bool}) = TruthTableConstraint(BitMatrix(tt))

# Get pin number from constraint
get_pin_num(c::TruthTableConstraint) = size(c.truth_table, 2)

# ============================================================================
# Unified Gadget Types
# ============================================================================

"""
    Gadget{M<:EnergyModel, T<:Real}

A gadget found by the search algorithm.

# Type Parameters
- `M`: Energy model type (RydbergModel or QUBOModel)
- `T`: Numeric type for weights

# Fields
- `constraint::GadgetConstraint`: The constraint this gadget satisfies
- `graph::SimpleGraph{Int}`: The graph structure
- `pins::Vector{Int}`: Pin vertex indices
- `vertex_weights::Vector{T}`: Vertex weights (hᵢ)
- `edge_weights::Vector{T}`: Edge weights (Jᵢⱼ), empty for RydbergModel
- `edge_list::Vector{Tuple{Int,Int}}`: Edge list corresponding to edge_weights
- `pos::Union{Nothing, Vector{Tuple{Float64, Float64}}}`: Vertex positions
"""
struct Gadget{M<:EnergyModel, T<:Real}
    constraint::GadgetConstraint
    graph::SimpleGraph{Int}
    pins::Vector{Int}
    vertex_weights::Vector{T}
    edge_weights::Vector{T}
    edge_list::Vector{Tuple{Int,Int}}
    pos::Union{Nothing, Vector{Tuple{Float64, Float64}}}
end

# Convenience constructors
function Gadget(::Type{RydbergModel}, constraint::GadgetConstraint, graph::SimpleGraph{Int}, 
                pins::Vector{Int}, vertex_weights::Vector{T}, 
                pos::Union{Nothing, Vector{Tuple{Float64, Float64}}}=nothing) where T<:Real
    Gadget{RydbergModel, T}(constraint, graph, pins, vertex_weights, T[], Tuple{Int,Int}[], pos)
end

function Gadget(::Type{QUBOModel}, constraint::GadgetConstraint, graph::SimpleGraph{Int}, 
                pins::Vector{Int}, vertex_weights::Vector{T}, edge_weights::Vector{T},
                edge_list::Vector{Tuple{Int,Int}},
                pos::Union{Nothing, Vector{Tuple{Float64, Float64}}}=nothing) where T<:Real
    Gadget{QUBOModel, T}(constraint, graph, pins, vertex_weights, edge_weights, edge_list, pos)
end

function Gadget(::Type{RydbergUnweightedModel}, constraint::GadgetConstraint, graph::SimpleGraph{Int}, 
                pins::Vector{Int}, vertex_weights::Vector{T}, 
                pos::Union{Nothing, Vector{Tuple{Float64, Float64}}}=nothing) where T<:Real
    Gadget{RydbergUnweightedModel, T}(constraint, graph, pins, vertex_weights, T[], Tuple{Int,Int}[], pos)
end

# Legacy compatibility: Gadget with truth table and weights only
function Gadget(ground_states::BitMatrix, graph::SimpleGraph{Int}, pins::Vector{Int}, 
                weights::Vector{T}, pos::Union{Nothing, Vector{Tuple{Float64, Float64}}}=nothing) where T<:Real
    Gadget{RydbergModel, T}(TruthTableConstraint(ground_states), graph, pins, weights, T[], Tuple{Int,Int}[], pos)
end

# Accessor for legacy compatibility
function Base.getproperty(g::Gadget, s::Symbol)
    if s == :weights
        return getfield(g, :vertex_weights)
    elseif s == :ground_states
        c = getfield(g, :constraint)
        return c.truth_table
    else
        return getfield(g, s)
    end
end

# ============================================================================
# State Space Functions
# ============================================================================

"""
    get_state_space(::Type{RydbergModel}, graph::SimpleGraph{Int}) -> Tuple{Vector{UInt32}, Int}

Get the state space for Rydberg model (Maximal Independent Sets).
"""
function get_state_space(::Type{RydbergModel}, graph::SimpleGraph{Int})
    return find_maximal_independent_sets(graph)
end

"""
    get_state_space(::Type{QUBOModel}, graph::SimpleGraph{Int}) -> Tuple{Vector{UInt32}, Int}

Get the state space for QUBO model (all 2^n states).
"""
function get_state_space(::Type{QUBOModel}, graph::SimpleGraph{Int})
    n = Graphs.nv(graph)
    n > MAX_SUPPORTED_VERTICES && error("Too many vertices: $n > $MAX_SUPPORTED_VERTICES")
    states = collect(UInt32(0):UInt32(2^n - 1))
    return (states, length(states))
end

"""
    get_state_space(::Type{RydbergUnweightedModel}, graph::SimpleGraph{Int})

Get the state space for unweighted Rydberg model (Maximal Independent Sets).
Same state space as RydbergModel.
"""
function get_state_space(::Type{RydbergUnweightedModel}, graph::SimpleGraph{Int})
    return find_maximal_independent_sets(graph)
end


"""
    _graph_hash(g::SimpleGraph{Int}) -> UInt64

Create a hash for graph structure for caching purposes.
"""
function _graph_hash(g::SimpleGraph{Int})::UInt64
    edges = sort!([(min(src(e), dst(e)), max(src(e), dst(e))) for e in Graphs.edges(g)])
    return hash((Graphs.nv(g), edges))
end

"""
    find_maximal_independent_sets(g)

Find all maximal independent sets of a graph using bit masks with caching.
"""
function find_maximal_independent_sets(g::SimpleGraph{Int})
    vertex_count = Graphs.nv(g)
    
    if vertex_count > MAX_SUPPORTED_VERTICES
        error("Graph has $(vertex_count) vertices, but maximum supported is $(MAX_SUPPORTED_VERTICES).")
    end
    
    graph_hash = _graph_hash(g)
    if haskey(MIS_CACHE, graph_hash)
        return MIS_CACHE[graph_hash]
    end
    
    complement_g = Graphs.complement(g)
    cliques = Graphs.maximal_cliques(complement_g)
    
    masks = Vector{UInt32}(undef, length(cliques))
    
    @inbounds for (i, clique) in enumerate(cliques)
        mask::UInt32 = 0
        for v in clique
            mask |= UInt32(1) << (v - 1)
        end
        masks[i] = mask
    end
    
    result = (masks, length(masks))
    MIS_CACHE[graph_hash] = result
    
    return result
end

# ============================================================================
# Constraint Matching Functions
# ============================================================================

"""
    match_constraint_to_states(states, constraint, pin_set)

Match constraint ground states to full graph states based on pin configuration.

# Returns
- `Vector{Vector{Int}}`: For each ground state pattern, indices of matching states in the state space
"""
function match_constraint_to_states(
    states::Vector{UInt32}, 
    constraint::TruthTableConstraint, 
    pin_set::Vector{Int}
)::Vector{Vector{Int}}
    
    tt = constraint.truth_table
    pin_num = size(tt, 2)
    length(pin_set) != pin_num && 
        error("Pin set length $(length(pin_set)) != constraint pin_num $(pin_num)")
    
    num_ground_states = size(tt, 1)
    result = Vector{Vector{Int}}(undef, num_ground_states)
    
    for gs_idx in 1:num_ground_states
        # Convert truth table row to target mask
        target_mask::UInt32 = 0
        for (bit_pos, val) in enumerate(tt[gs_idx, :])
            if val
                target_mask |= UInt32(1) << (bit_pos - 1)
            end
        end
        
        matches = Int[]
        for (state_idx, state) in enumerate(states)
            extracted::UInt32 = 0
            for (bit_pos, pin) in enumerate(pin_set)
                extracted |= ((state >> (pin - 1)) & 0x1) << (bit_pos - 1)
            end
            
            if extracted == target_mask
                push!(matches, state_idx)
            end
        end
        result[gs_idx] = matches
    end
    
    return result
end

# Legacy compatibility alias
match_rows_by_pinset(masks::Vector{UInt32}, truth_table::BitMatrix, pin_set::Vector{Int}) = 
    match_constraint_to_states(masks, TruthTableConstraint(truth_table), pin_set)

# ============================================================================
# Weight Finding Functions
# ============================================================================

"""
    solve_weights(model_type, states, target_indices_all, graph, pin_set, optimizer; kwargs...)

Unified weight solver that works for both Rydberg and QUBO models.
"""
function solve_weights(
    ::Type{M},
    states::Vector{UInt32},
    target_indices_all::Vector{Vector{Int}},
    graph::SimpleGraph{Int},
    pin_set::Vector{Int},
    optimizer,
    env=nothing,
    objective=nothing,
    allow_defect::Bool=false,
    max_samples::Int=1000,
    check_connectivity::Bool=true
) where M <: EnergyModel
    if optimizer === nothing && M !== RydbergUnweightedModel
        error("Optimizer must be provided for $(M). Only RydbergUnweightedModel can run without an optimizer.")
    end

    any(isempty, target_indices_all) && return nothing
    
    vertex_num = Graphs.nv(graph)
    edge_list = [(src(e), dst(e)) for e in Graphs.edges(graph)]
    all_state_indices = collect(1:length(states))
    
    total_combinations = prod(length, target_indices_all)
    if total_combinations > max_samples
        @warn "Too many combinations ($total_combinations), sampling first $max_samples"
        return _solve_weights_sampled(M, states, target_indices_all, graph, vertex_num, edge_list,
                                      pin_set, optimizer, env, objective, allow_defect, 
                                      max_samples, all_state_indices, check_connectivity)
    end
    
    target_indices_set = Vector{Int}(undef, length(target_indices_all))
    
    sample_count = 0
    for target_indices in Iterators.product(target_indices_all...)
        sample_count += 1
        if sample_count > max_samples
            break
        end
        
        copyto!(target_indices_set, target_indices)
        wrong_indices = setdiff(all_state_indices, target_indices_set)
        
        target_states = [states[i] for i in target_indices_set]
        wrong_states = [states[i] for i in wrong_indices]

        result = _find_weights(M, vertex_num, edge_list, pin_set, target_states, wrong_states, 
                               optimizer, env, objective, allow_defect, graph, check_connectivity)
        if result !== nothing
            return result
        end
    end

    return nothing
end

"""
    _solve_weights_sampled(...)

Optimized version for large combination spaces using random sampling.
"""
function _solve_weights_sampled(
    ::Type{M},
    states::Vector{UInt32},
    target_indices_all::Vector{Vector{Int}},
    graph::SimpleGraph{Int},
    vertex_num::Int,
    edge_list::Vector{Tuple{Int,Int}},
    pin_set::Vector{Int},
    optimizer,
    env,
    objective,
    allow_defect::Bool,
    max_samples::Int,
    all_state_indices::Vector{Int},
    check_connectivity::Bool=true
) where M <: EnergyModel
    target_indices_set = Vector{Int}(undef, length(target_indices_all))
    
    for _ in 1:max_samples
        for (i, indices) in enumerate(target_indices_all)
            target_indices_set[i] = rand(indices)
        end
        
        wrong_indices = setdiff(all_state_indices, target_indices_set)
        
        target_states = [states[i] for i in target_indices_set]
        wrong_states = [states[i] for i in wrong_indices]

        result = _find_weights(M, vertex_num, edge_list, pin_set, target_states, wrong_states, 
                               optimizer, env, objective, allow_defect, graph, check_connectivity)
        if result !== nothing
            return result
        end
    end
    
    return nothing
end

"""
    check_connectivity_after_removal(graph, vertices_to_remove)

Check if a graph remains connected after removing specified vertices.
"""
function check_connectivity_after_removal(graph::SimpleGraph{Int}, vertices_to_remove::Vector{Int})
    if isempty(vertices_to_remove)
        return Graphs.is_connected(graph)
    end
    
    vertices_to_remove_1based = vertices_to_remove .+ 1
    all_vertices = collect(1:Graphs.nv(graph))
    remaining_vertices = setdiff(all_vertices, vertices_to_remove_1based)
    
    if isempty(remaining_vertices)
        return false
    end
    
    subgraph, _ = Graphs.induced_subgraph(graph, remaining_vertices)
    return Graphs.is_connected(subgraph)
end

"""
    _find_weights(::Type{RydbergModel}, ...) -> Union{Nothing, Vector{Float64}}

Find vertex weights for Rydberg model.
"""
function _find_weights(
    ::Type{RydbergModel},
    vertex_num::Int, 
    edge_list::Vector{Tuple{Int,Int}}, 
    pin_set::Vector{Int}, 
    target_states::Vector{UInt32}, 
    wrong_states::Vector{UInt32}, 
    optimizer, 
    env, 
    objective, 
    allow_defect::Bool, 
    graph::SimpleGraph{Int},
    check_connectivity::Bool=true
)
    opt = isnothing(env) ? optimizer() : optimizer(env)
    model = direct_model(opt)
    set_silent(model)
    set_string_names_on_creation(model, false)

    if allow_defect
        x = @variable(model, [i=1:vertex_num], lower_bound = i in pin_set ? 1 : 0, Int)
    else
        @variable(model, x[1:vertex_num] >= 1, Int)
    end

    @variable(model, C)

    # Target states must have equal energy
    for m in target_states
        @constraint(model, sum(((m >> (v - 1)) & 0x1) * x[v] for v in 1:vertex_num) == C)
    end

    # Wrong states must have higher energy (for MIS, we want lower energy for ground states)
    for m in wrong_states
        @constraint(model, sum(((m >> (v - 1)) & 0x1) * x[v] for v in 1:vertex_num) <= C - DEFAULT_EPSILON)
    end

    if objective !== nothing
        @objective(model, Min, objective(x))
    end

    optimize!(model)
    if is_solved_and_feasible(model)
        weights = [value(x[v]) for v in 1:vertex_num]
        
        if check_connectivity
            zero_weight_vertices = findall(w -> abs(w) < 1e-6, weights)
            
            if !isempty(zero_weight_vertices)
                zero_vertices_0based = zero_weight_vertices .- 1
                
                if !check_connectivity_after_removal(graph, zero_vertices_0based)
                    @info "Solution found but graph becomes disconnected after removing zero-weight vertices, rejecting"
                    return nothing
                end
            end
        end
        
        @info "found a valid Rydberg solution"
        @show weights
        return weights
    end
    return nothing
end

"""
    _find_weights(::Type{QUBOModel}, ...) -> Union{Nothing, Tuple{Vector{Float64}, Vector{Float64}}}

Find vertex and edge weights for QUBO model.
"""
function _find_weights(
    ::Type{QUBOModel},
    vertex_num::Int, 
    edge_list::Vector{Tuple{Int,Int}}, 
    pin_set::Vector{Int}, 
    target_states::Vector{UInt32}, 
    wrong_states::Vector{UInt32}, 
    optimizer, 
    env, 
    objective, 
    allow_defect::Bool, 
    graph::SimpleGraph{Int},
    check_connectivity::Bool=true
)
    opt = isnothing(env) ? optimizer() : optimizer(env)
    model = direct_model(opt)
    set_silent(model)
    set_string_names_on_creation(model, false)

    # Vertex weights (h_i)
    @variable(model, h[1:vertex_num], Int)
    
    # Edge weights (J_ij)
    edge_num = length(edge_list)
    @variable(model, J[1:edge_num], Int)

    @variable(model, C)

    # Helper function to compute energy of a state
    function state_energy(state::UInt32)
        vertex_energy = sum(((state >> (v - 1)) & 0x1) * h[v] for v in 1:vertex_num)
        edge_energy = sum(
            ((state >> (edge_list[e][1] - 1)) & 0x1) * 
            ((state >> (edge_list[e][2] - 1)) & 0x1) * J[e] 
            for e in 1:edge_num;
            init=0
        )
        return vertex_energy + edge_energy
    end

    # Target states must have equal energy C
    for state in target_states
        @constraint(model, state_energy(state) == C)
    end

    # Wrong states must have higher energy
    for state in wrong_states
        @constraint(model, state_energy(state) <= C - DEFAULT_EPSILON)
    end

    if objective !== nothing
        @objective(model, Min, objective(h, J))
    end

    optimize!(model)
    if is_solved_and_feasible(model)
        vertex_weights = [value(h[v]) for v in 1:vertex_num]
        edge_weights = edge_num > 0 ? [value(J[e]) for e in 1:edge_num] : Float64[]
        
        for i in pin_set
            if vertex_weights[i] == 0
                return nothing
            end
        end
        if check_connectivity
            zero_weight_vertices = findall(w -> abs(w) < 1e-6, vertex_weights)
            
            if !isempty(zero_weight_vertices)
                zero_vertices_0based = zero_weight_vertices .- 1
                
                if !check_connectivity_after_removal(graph, zero_vertices_0based)
                    @info "QUBO solution found but graph becomes disconnected, rejecting"
                    return nothing
                end
            end
        end
        
        @info "found a valid QUBO solution"
        @show vertex_weights, edge_weights
        return (vertex_weights, edge_weights)
    end
    return nothing
end

"""
    _find_weights(::Type{RydbergUnweightedModel}, ...) -> Union{Nothing, Vector{Float64}}

Verify gadget using reduced α-tensor for the unweighted Rydberg model.

Computes the reduced α-tensor α[s] for each boundary (pin) configuration s,
then determines ground states as the configurations maximizing:
    total(s) = α[s] + count_ones(s)

Checks that these ground configurations match exactly the pin configurations
of the target states (which correspond to the truth table).

No optimizer is needed — this is a pure combinatorial verification.

Reference: Liu et al., "Computer-assisted gadget design and problem reduction
of unweighted maximum independent set", PRX Quantum 4, 010316 (2023)
"""
function _find_weights(
    ::Type{RydbergUnweightedModel},
    vertex_num::Int, 
    edge_list::Vector{Tuple{Int,Int}}, 
    pin_set::Vector{Int}, 
    target_states::Vector{UInt32}, 
    wrong_states::Vector{UInt32}, 
    optimizer,  # Not used
    env, 
    objective, 
    allow_defect::Bool, 
    graph::SimpleGraph{Int},
    check_connectivity::Bool=true
)
    # Use α-tensor verification
    result = verify_gadget_via_alpha_tensor(graph, pin_set, target_states)
    
    if result === nothing
        return nothing
    end
    
    weights, max_total = result
    @info "found a valid RydbergUnweighted solution via α-tensor (MIS size = $max_total)"
    return weights
end

# ============================================================================
# Unified Search Interface
# ============================================================================

"""
    make_filter(model_type, constraint, optimizer, env; kwargs...)

Create a filter closure for graph search that checks if a graph can implement the given constraint.

# Arguments
- `model_type::Type{<:EnergyModel}`: RydbergModel or QUBOModel
- `constraint::GadgetConstraint`: The target constraint to implement
- `optimizer`: Optimization solver for weight finding
- `env`: Environment for the optimizer

# Keywords
- `connected::Bool=false`: Whether to require connected graphs only
- `objective=nothing`: Objective function for optimization
- `allow_defect::Bool=false`: Whether to allow zero vertex weights
- `max_samples::Int=1000`: Maximum samples for enumeration
- `pin_candidates::Union{Nothing, Vector{Vector{Int}}}=nothing`: Candidate pin combinations
- `check_connectivity::Bool=true`: Whether to check graph connectivity

# Returns
- `Function`: Filter function that takes (graph, pos, pin_set) and returns Gadget or nothing
"""
function make_filter(
    ::Type{M},
    constraint::GadgetConstraint, 
    optimizer, 
    env; 
    connected::Bool=false, 
    objective=nothing, 
    allow_defect::Bool=false, 
    max_samples::Int=1000, 
    pin_candidates::Union{Nothing, Vector{Vector{Int}}}=nothing, 
    check_connectivity::Bool=true
) where M <: EnergyModel
  
    return function(graph::SimpleGraph{Int}, pos::Vector{Tuple{Float64, Float64}}, pin_set::Union{Nothing, Vector{Int}}=nothing)
        if !connected
            Graphs.is_connected(graph) || return nothing
        end

        vertex_num = Graphs.nv(graph)
        states, _ = get_state_space(M, graph)
        edge_list = [(src(e), dst(e)) for e in Graphs.edges(graph)]

        if pin_set === nothing
            pin_set = collect(1:vertex_num)
        end
        
        pin_num = get_pin_num(constraint)
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
            all_candidates = collect(Combinatorics.combinations(pin_set, pin_num))
        end
        
        for candidate in all_candidates
            target_indices_all = match_constraint_to_states(states, constraint, candidate)
            result = solve_weights(M, states, target_indices_all, graph, candidate, 
                                   optimizer, env, objective, allow_defect, max_samples, check_connectivity)
            if result !== nothing
                if M === QUBOModel
                    vertex_weights, edge_weights = result
                    return Gadget(QUBOModel, constraint, graph, candidate, vertex_weights, edge_weights, edge_list, pos)
                else  # RydbergModel or RydbergUnweightedModel
                    return Gadget(M, constraint, graph, candidate, result, pos)
                end
            end
        end
        return nothing
    end
end

# Legacy compatibility: make_filter with truth_table
function make_filter(truth_table::BitMatrix, optimizer, env; kwargs...)
    make_filter(RydbergModel, TruthTableConstraint(truth_table), optimizer, env; kwargs...)
end

"""
    search_gadgets(model_type, loader, constraints; kwargs...)

Unified search function for both Rydberg and QUBO gadgets.

# Arguments
- `model_type::Type{<:EnergyModel}`: RydbergModel or QUBOModel
- `loader::GraphLoader`: The graph loader containing candidate graphs
- `constraints::Vector{<:GadgetConstraint}`: Constraints to search for

# Keywords
- `optimizer`: Optimization solver to use for weight finding
- `env=nothing`: Environment for the optimizer
- `connected::Bool=false`: Whether to require connected graphs only
- `objective=nothing`: Objective function for optimization
- `allow_defect::Bool=false`: Whether to allow zero vertex weights
- `limit=nothing`: Maximum number of graphs to search
- `max_samples::Int=100`: Maximum samples for weight enumeration
- `save_path::String="results.json"`: Path to save intermediate results
- `pin_candidates::Union{Nothing, Vector{Vector{Int}}}=nothing`: Candidate pin combinations
- `check_connectivity::Bool=true`: Whether to check graph connectivity
- `max_result_num::Int=1`: Maximum number of results per constraint

# Returns
- `Tuple{Vector{Vector{Gadget}}, Vector{<:GadgetConstraint}}`: Found gadgets grouped by constraint and failed constraints
"""
function search_gadgets(
    ::Type{M},
    loader::GraphLoader,
    constraints::Vector{C};
    optimizer=nothing,
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
) where {M <: EnergyModel, C <: GadgetConstraint}
    
    results = Vector{Vector{Gadget}}(undef, length(constraints))
    for i in 1:length(constraints)
        results[i] = Gadget[]
    end
    failed_constraints = C[]
    total_found = 0

    start_time = time()
    initial_cache_size = length(MIS_CACHE)

    for (i, constraint) in enumerate(constraints)
        constraint_start = time()
        model_name = M === RydbergModel ? "Rydberg" : M === RydbergUnweightedModel ? "RydbergUnweighted (α-tensor)" : "QUBO"
        @info "[$model_name] Searching for constraint $(i-1) [limit=$max_result_num]"
        
        filter_fn = make_filter(M, constraint, optimizer, env;
                                connected=connected,
                                objective=objective,
                                allow_defect=allow_defect,
                                max_samples=max_samples,
                                pin_candidates=pin_candidates,
                                check_connectivity=check_connectivity)
        
        gadgets = find_matching_gadget(loader; filter=filter_fn, limit=limit, max_results=max_result_num)
        
        constraint_time = time() - constraint_start
        @info "Constraint $(i-1) processed in $(round(constraint_time, digits=2))s, found $(length(gadgets)) gadgets"
        
        if !isempty(gadgets)
            results[i] = gadgets
            total_found += length(gadgets)
            
            all_gadgets = vcat(results...)
            save_results_to_json(all_gadgets, save_path)
        else
            push!(failed_constraints, constraint)
        end
        
        if length(MIS_CACHE) > initial_cache_size + 1000
            @info "Cache growing large ($(length(MIS_CACHE)) entries)"
        end
    end

    total_time = time() - start_time
    cache_hits = length(MIS_CACHE) - initial_cache_size
    @info "Search completed in $(round(total_time, digits=2))s. Cache gained $cache_hits entries."

    return results, failed_constraints
end

# ============================================================================
# Legacy Compatibility Functions
# ============================================================================

"""
    search_by_truth_tables(loader, truth_tables; kwargs...)

Search for Rydberg gadgets by truth tables (legacy interface).
"""
function search_by_truth_tables(
    loader::GraphLoader,
    truth_tables::Vector{BitMatrix};
    kwargs...
)
    constraints = [TruthTableConstraint(tt) for tt in truth_tables]
    results, failed = search_gadgets(RydbergModel, loader, constraints; kwargs...)
    failed_tt = [f.truth_table for f in failed]
    return results, failed_tt
end

"""
    find_matching_gadget(loader; filter=nothing, limit=nothing, keys_range=nothing, max_results=nothing)

Find gadgets matching a given filter function from the graph loader.
"""
function find_matching_gadget(loader::GraphLoader; filter=nothing, limit::Union{Int,Nothing}=nothing, keys_range::Union{Nothing, Vector{Int}}=nothing, max_results::Union{Int,Nothing}=nothing)
    keys_raw = keys_range === nothing ? collect(keys(loader)) : keys_range
    keys_to_search = isa(keys_raw[1], Int) ? keys_raw : parse.(Int, keys_raw)
    total = limit === nothing ? length(keys_to_search) : min(length(keys_to_search), limit)

    results = Gadget[]

    @showprogress for key in Iterators.take(keys_to_search, total)
        g = loader[key]
        result = filter(g, loader.layout[key], loader.pinset)
        if result !== nothing
            push!(results, result)
            if max_results !== nothing && length(results) >= max_results
                @info "Early termination: found $(length(results)) results (max_results=$max_results)"
                break
            end
        end
    end
    return results
end

# Legacy solve_weight_enumerate for compatibility
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
    if graph === nothing
        error("Graph must be provided for solve_weight_enumerate")
    end
    result = solve_weights(RydbergModel, mis_result, target_mis_indices_all, graph, pin_set,
                           optimizer, env, objective, allow_defect, max_samples, check_connectivity)
    return result === nothing ? Float64[] : result
end
