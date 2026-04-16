"""
    RydbergModel <: EnergyModel

Energy model for Rydberg atom systems.
- State space: Maximal Independent Sets (MIS)
- Energy: E(σ) = Σᵢ hᵢσᵢ (vertex weights only)
"""
struct RydbergModel <: EnergyModel end

# ============================================================================
# State Space
# ============================================================================

"""
    get_state_space(::Type{RydbergModel}, graph) -> Tuple{Vector{UInt32}, Int}

Get the state space for Rydberg model (Maximal Independent Sets).
"""
function get_state_space(::Type{RydbergModel}, graph::SimpleGraph{Int})
    return find_maximal_independent_sets(graph)
end

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
# Weight Finding
# ============================================================================

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
    model = _make_model(optimizer, env)
    if allow_defect
        x = @variable(model, [i=1:vertex_num], lower_bound = i in pin_set ? 1 : 0, Int)
    else
        @variable(model, x[1:vertex_num] >= 1, Int)
    end
    @variable(model, C)
    for m in target_states
        @constraint(model, sum(((m >> (v - 1)) & 0x1) * x[v] for v in 1:vertex_num) == C)
    end
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

# ============================================================================
# Rydberg Gadget Constructor
# ============================================================================

function Gadget(::Type{RydbergModel}, constraint::GadgetConstraint, graph::SimpleGraph{Int},
                pins::Vector{Int}, vertex_weights::Vector{T},
                pos::Union{Nothing, Vector{Tuple{Float64, Float64}}}=nothing) where T<:Real
    Gadget{RydbergModel, T}(constraint, graph, pins, vertex_weights, T[], Tuple{Int,Int}[], pos)
end

function _build_gadget(::Type{RydbergModel}, constraint, graph, candidate, result, edge_list, pos)
    return Gadget(RydbergModel, constraint, graph, candidate, result, pos)
end

# Legacy constructor
function Gadget(ground_states::BitMatrix, graph::SimpleGraph{Int}, pins::Vector{Int},
                weights::Vector{T}, pos::Union{Nothing, Vector{Tuple{Float64, Float64}}}=nothing) where T<:Real
    Gadget{RydbergModel, T}(TruthTableConstraint(ground_states), graph, pins, weights, T[], Tuple{Int,Int}[], pos)
end

# ============================================================================
# Legacy Compatibility
# ============================================================================

function make_filter(truth_table::BitMatrix, optimizer, env; kwargs...)
    make_filter(RydbergModel, TruthTableConstraint(truth_table), optimizer, env; kwargs...)
end

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

match_rows_by_pinset(masks::Vector{UInt32}, truth_table::BitMatrix, pin_set::Vector{Int}) =
    match_constraint_to_states(masks, TruthTableConstraint(truth_table), pin_set)

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
