"""
    QUBOModel <: EnergyModel

Energy model for general QUBO (Quadratic Unconstrained Binary Optimization).
- State space: All 2^n binary configurations
- Energy: E(σ) = Σᵢ hᵢσᵢ + Σ₍ᵢ,ⱼ₎ Jᵢⱼσᵢσⱼ (vertex + edge weights)
"""
struct QUBOModel <: EnergyModel end

# ============================================================================
# State Space
# ============================================================================

"""
    get_state_space(::Type{QUBOModel}, graph) -> Tuple{Vector{UInt32}, Int}

Get the state space for QUBO model (all 2^n states).
"""
function get_state_space(::Type{QUBOModel}, graph::SimpleGraph{Int})
    n = Graphs.nv(graph)
    n > MAX_SUPPORTED_VERTICES && error("Too many vertices: $n > $MAX_SUPPORTED_VERTICES")
    states = collect(UInt32(0):UInt32(2^n - 1))
    return (states, length(states))
end

# ============================================================================
# Weight Finding
# ============================================================================

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
    model = _make_model(optimizer, env)
    W_MAX = 10
    @variable(model, -W_MAX <= h[1:vertex_num] <= W_MAX, Int)
    edge_num = length(edge_list)
    @variable(model, -W_MAX <= J[1:edge_num] <= W_MAX, Int)
    @variable(model, z_edge[1:edge_num], Bin)
    for e in 1:edge_num
        @constraint(model, J[e] >= 1 - 2 * W_MAX * z_edge[e])
        @constraint(model, J[e] <= -1 + 2 * W_MAX * (1 - z_edge[e]))
    end
    @variable(model, C)
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
    for state in target_states
        @constraint(model, state_energy(state) == C)
    end
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
        @info "found a valid QUBO solution"
        @show vertex_weights, edge_weights
        return (vertex_weights, edge_weights)
    end
    return nothing
end

# ============================================================================
# QUBO Gadget Constructor
# ============================================================================

function Gadget(::Type{QUBOModel}, constraint::GadgetConstraint, graph::SimpleGraph{Int},
                pins::Vector{Int}, vertex_weights::Vector{T}, edge_weights::Vector{T},
                edge_list::Vector{Tuple{Int,Int}},
                pos::Union{Nothing, Vector{Tuple{Float64, Float64}}}=nothing) where T<:Real
    Gadget{QUBOModel, T}(constraint, graph, pins, vertex_weights, edge_weights, edge_list, pos)
end

function _build_gadget(::Type{QUBOModel}, constraint, graph, candidate, result, edge_list, pos)
    vertex_weights, edge_weights = result
    return Gadget(QUBOModel, constraint, graph, candidate, vertex_weights, edge_weights, edge_list, pos)
end
