# α-Tensor computation for gadget verification
# Based on the reduced α-tensor framework from:
#   Liu et al., "Computer-assisted gadget design and problem reduction 
#   of unweighted maximum independent set", PRX Quantum 4, 010316 (2023)
#
# The α-tensor α[s] gives the maximum independent set size in the interior
# of a graph, given a fixed boundary (pin) configuration s.
# Ground states are the boundary configurations that maximize:
#   total(s) = α[s] + count_ones(s)

"""Sentinel value for infeasible boundary configurations (e.g., adjacent pins both selected)."""
const INFEASIBLE_ALPHA = -1

"""
    compute_reduced_alpha_tensor(graph::SimpleGraph{Int}, pins::Vector{Int})

Compute the reduced α-tensor for a graph with designated pin (boundary) vertices.

For each boundary configuration (which pins are in the independent set), 
compute the maximum cardinality of an independent set in the interior vertices,
subject to the constraint that selected pins' neighbors in the interior are forbidden.

Infeasible configurations (where selected pins are adjacent) get value `INFEASIBLE_ALPHA`.

# Arguments
- `graph`: The graph structure
- `pins`: Indices of boundary/pin vertices

# Returns
- `Dict{UInt32, Int}`: Mapping from boundary configuration to max interior MIS size
"""
function compute_reduced_alpha_tensor(
    graph::SimpleGraph{Int},
    pins::Vector{Int}
)
    n_pins = length(pins)
    n_vertices = Graphs.nv(graph)
    
    # Get interior vertices (non-pins)
    all_vertices = collect(1:n_vertices)
    interior = setdiff(all_vertices, pins)
    
    alpha = Dict{UInt32, Int}()
    
    # Enumerate all 2^n_pins boundary configurations
    for boundary_mask in 0:(2^n_pins - 1)
        boundary_config = UInt32(boundary_mask)
        
        # Find max MIS in interior with this boundary fixed
        max_interior_size = compute_max_interior_mis(
            graph, interior, pins, boundary_config
        )
        
        alpha[boundary_config] = max_interior_size
    end
    
    return alpha
end

"""
    compute_maximum_independent_set_size(graph::SimpleGraph{Int})

Compute the size of the maximum independent set (not just maximal) in a graph.

This uses exhaustive enumeration for small graphs. For the interior vertices
in gadget verification, graphs are typically small enough for this approach.

# Returns
- Maximum size of an independent set
"""
function compute_maximum_independent_set_size(graph::SimpleGraph{Int})
    n = Graphs.nv(graph)
    if n == 0
        return 0
    end
    
    # For small graphs, use exhaustive enumeration
    max_size = 0
    
    # Enumerate all possible subsets
    for mask in 0:(2^n - 1)
        # Extract vertices in this subset
        vertices_in_set = Int[]
        for v in 1:n
            if ((mask >> (v - 1)) & 0x1) == 1
                push!(vertices_in_set, v)
            end
        end
        
        # Check if this subset is an independent set
        is_independent = true
        for i in 1:length(vertices_in_set)
            for j in (i+1):length(vertices_in_set)
                if Graphs.has_edge(graph, vertices_in_set[i], vertices_in_set[j])
                    is_independent = false
                    break
                end
            end
            !is_independent && break
        end
        
        if is_independent
            max_size = max(max_size, length(vertices_in_set))
        end
    end
    
    return max_size
end

"""
    compute_max_interior_mis(graph, interior, pins, boundary_config)

Compute the maximum size of an independent set in interior vertices,
given that boundary vertices are fixed according to boundary_config.

Returns `INFEASIBLE_ALPHA` if the boundary configuration is infeasible
(i.e., two adjacent pins are both selected).

# Arguments
- `graph`: The graph
- `interior`: Interior vertex indices
- `pins`: Pin vertex indices  
- `boundary_config`: Bit mask indicating which pins are in the set
  (bit i-1 corresponds to pins[i])

# Returns
- Maximum cardinality of interior independent set compatible with boundary,
  or `INFEASIBLE_ALPHA` if boundary is infeasible
"""
function compute_max_interior_mis(
    graph::SimpleGraph{Int},
    interior::Vector{Int},
    pins::Vector{Int},
    boundary_config::UInt32
)
    # Determine which pins are selected
    selected_pins = Int[]
    for (i, pin) in enumerate(pins)
        if ((boundary_config >> (i - 1)) & 0x1) == 1
            push!(selected_pins, pin)
        end
    end
    
    # Check boundary feasibility: selected pins must form an independent set
    for i in 1:length(selected_pins)
        for j in (i+1):length(selected_pins)
            if Graphs.has_edge(graph, selected_pins[i], selected_pins[j])
                return INFEASIBLE_ALPHA  # Adjacent pins both selected → infeasible
            end
        end
    end
    
    # Find vertices in interior that are neighbors of selected pins
    # These vertices are forbidden in the interior MIS
    forbidden = Set{Int}()
    for pin in selected_pins
        for neighbor in Graphs.neighbors(graph, pin)
            if neighbor ∈ interior
                push!(forbidden, neighbor)
            end
        end
    end
    
    # Create induced subgraph on feasible interior vertices
    feasible_interior = setdiff(interior, forbidden)
    
    if isempty(feasible_interior)
        return 0
    end
    
    # Compute maximum independent set size on the induced subgraph
    subgraph, vmap = Graphs.induced_subgraph(graph, feasible_interior)
    return compute_maximum_independent_set_size(subgraph)
end

"""
    find_ground_configs(α::Dict{UInt32, Int}, n_pins::Int)

From an α-tensor, find the boundary configurations that are ground states.

Ground states maximize total(s) = α[s] + count_ones(s) among all feasible configs.

# Returns
- `(ground_configs::Set{UInt32}, max_total::Int)`: The set of ground state
  boundary configurations and the maximum total IS size.
  Returns `(empty set, -1)` if no feasible configuration exists.
"""
function find_ground_configs(α::Dict{UInt32, Int}, n_pins::Int)
    max_total = typemin(Int)
    
    for config in UInt32(0):UInt32(2^n_pins - 1)
        α_val = get(α, config, INFEASIBLE_ALPHA)
        α_val == INFEASIBLE_ALPHA && continue
        total = α_val + count_ones(config)
        max_total = max(max_total, total)
    end
    
    if max_total == typemin(Int)
        return Set{UInt32}(), -1
    end
    
    ground_configs = Set{UInt32}()
    for config in UInt32(0):UInt32(2^n_pins - 1)
        α_val = get(α, config, INFEASIBLE_ALPHA)
        α_val == INFEASIBLE_ALPHA && continue
        if α_val + count_ones(config) == max_total
            push!(ground_configs, config)
        end
    end
    
    return ground_configs, max_total
end

"""
    extract_pin_config(state::UInt32, pins::Vector{Int})

Extract the configuration of pin vertices from a full state.

# Arguments
- `state`: Full vertex configuration as bit mask
- `pins`: Pin vertex indices (1-based)

# Returns
- `UInt32`: Pin configuration as bit mask (pins[i] → bit i-1)
"""
function extract_pin_config(state::UInt32, pins::Vector{Int})
    config = UInt32(0)
    for (i, pin) in enumerate(pins)
        if ((state >> (pin - 1)) & 0x1) == 1
            config |= UInt32(1) << (i - 1)
        end
    end
    return config
end

"""
    check_alpha_equivalence(α1::Dict{UInt32, Int}, α2::Dict{UInt32, Int})

Check if two α-tensors are equivalent up to a constant.

Two α-tensors are equivalent if there exists a constant c such that:
    α2[config] = α1[config] + c  for all feasible configurations

Infeasible configurations (value `INFEASIBLE_ALPHA`) are skipped.

# Returns
- `(is_equivalent::Bool, constant::Union{Int, Nothing})`
"""
function check_alpha_equivalence(
    α1::Dict{UInt32, Int},
    α2::Dict{UInt32, Int}
)
    # Must have same keys (same boundary configurations)
    if Set(keys(α1)) != Set(keys(α2))
        return false, nothing
    end
    
    if isempty(α1)
        return true, 0
    end
    
    # Find the first feasible config to compute c
    c = nothing
    for config in keys(α1)
        v1 = α1[config]
        v2 = α2[config]
        # Skip infeasible configs
        if v1 == INFEASIBLE_ALPHA || v2 == INFEASIBLE_ALPHA
            # Both must be infeasible or equivalence fails
            if v1 != v2
                # One feasible, one not → not equivalent
                if (v1 == INFEASIBLE_ALPHA) != (v2 == INFEASIBLE_ALPHA)
                    return false, nothing
                end
            end
            continue
        end
        c = v2 - v1
        break
    end
    
    if c === nothing
        # All configs are infeasible
        return true, 0
    end
    
    # Verify constant holds for all feasible configurations
    for config in keys(α1)
        v1 = α1[config]
        v2 = α2[config]
        if v1 == INFEASIBLE_ALPHA && v2 == INFEASIBLE_ALPHA
            continue  # Both infeasible → OK
        elseif v1 == INFEASIBLE_ALPHA || v2 == INFEASIBLE_ALPHA
            return false, nothing  # One feasible, one not
        end
        if v2 - v1 != c
            return false, nothing
        end
    end
    
    return true, c
end

"""
    verify_gadget_via_alpha_tensor(
        graph::SimpleGraph{Int},
        pins::Vector{Int},
        target_states::Vector{UInt32}
    ) -> Union{Nothing, Tuple{Vector{Float64}, Int}}

Verify if a graph can serve as a gadget using α-tensor ground state analysis.

Computes the reduced α-tensor, finds which boundary (pin) configurations
maximize total(s) = α[s] + count_ones(s), and checks if these "ground configs"
match exactly the pin configurations of the target states.

# Returns
- `nothing` if gadget is invalid
- `(weights::Vector{Float64}, max_total::Int)` if valid,
  where weights are uniform (all 1.0) and max_total is the MIS size
"""
function verify_gadget_via_alpha_tensor(
    graph::SimpleGraph{Int},
    pins::Vector{Int},
    target_states::Vector{UInt32}
)
    # 1. Compute α-tensor
    α = compute_reduced_alpha_tensor(graph, pins)
    n_pins = length(pins)
    
    # 2. Find ground configs from α-tensor
    ground_configs, max_total = find_ground_configs(α, n_pins)
    
    if max_total < 0
        return nothing  # No feasible config
    end
    
    # 3. Extract target pin configs
    target_pin_configs = Set{UInt32}()
    for state in target_states
        push!(target_pin_configs, extract_pin_config(state, pins))
    end
    
    # 4. Ground configs must match exactly
    if ground_configs != target_pin_configs
        return nothing
    end
    
    # Valid gadget - return uniform weights
    n_vertices = Graphs.nv(graph)
    weights = ones(Float64, n_vertices)
    
    return (weights, max_total)
end
