# α-Tensor computation for gadget verification
# Based on UnitDiskMapping.jl's reduced α-tensor framework
#
# Reference: Liu et al., "Computer-assisted gadget design and problem reduction 
# of unweighted maximum independent set", PRX Quantum 4, 010316 (2023)

"""
    compute_reduced_alpha_tensor(graph::SimpleGraph{Int}, pins::Vector{Int})

Compute the reduced α-tensor for a graph with designated pin (boundary) vertices.

For each boundary configuration (which pins are in the independent set), 
compute the maximum cardinality of an independent set in the interior vertices.

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
        # Construct the fixed boundary configuration
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
    compute_max_interior_mis(graph, interior, pins, boundary_config)

Compute the maximum size of an independent set in interior vertices,
given that boundary vertices are fixed according to boundary_config.

# Arguments
- `graph`: The graph
- `interior`: Interior vertex indices
- `pins`: Pin vertex indices  
- `boundary_config`: Bit mask indicating which pins are in the set

# Returns
- Maximum cardinality of interior independent set compatible with boundary
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
    
    # Compute MIS on the induced subgraph
    subgraph, vmap = Graphs.induced_subgraph(graph, feasible_interior)
    
    # Find all maximal independent sets of the subgraph
    mis_masks, _ = find_maximal_independent_sets(subgraph)
    
    # Return the maximum cardinality
    return maximum(count_ones, mis_masks; init=0)
end

"""
    check_alpha_equivalence(α1::Dict{UInt32, Int}, α2::Dict{UInt32, Int})

Check if two α-tensors are equivalent up to a constant.

Two α-tensors are equivalent if there exists a constant c such that:
    α2[config] = α1[config] + c  for all configurations

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
    
    # Compute the constant from first configuration
    first_config = first(keys(α1))
    c = α2[first_config] - α1[first_config]
    
    # Verify constant holds for all configurations
    for config in keys(α1)
        if α2[config] - α1[config] != c
            return false, nothing
        end
    end
    
    return true, c
end

"""
    infer_pattern_alpha(target_states::Vector{UInt32}, pins::Vector{Int})

Infer the pattern's α-tensor from the target MIS states.

Given the target states that should be ground states, reconstruct what
the pattern's reduced α-tensor must be.

# Arguments
- `target_states`: MIS configurations that should be ground states
- `pins`: Pin vertex indices

# Returns
- `Dict{UInt32, Int}`: Inferred α-tensor for the pattern
"""
function infer_pattern_alpha(
    target_states::Vector{UInt32},
    pins::Vector{Int}
)
    alpha = Dict{UInt32, Int}()
    
    # Group target states by their pin configuration
    for state in target_states
        # Extract pin configuration from state
        pin_config = extract_pin_config(state, pins)
        
        # Count total vertices in this state (cardinality)
        total_size = count_ones(state)
        
        # Interior size = total - number of selected pins
        n_selected_pins = count_ones(pin_config)
        interior_size = total_size - n_selected_pins
        
        # Update alpha with maximum seen for this pin config
        if haskey(alpha, pin_config)
            alpha[pin_config] = max(alpha[pin_config], interior_size)
        else
            alpha[pin_config] = interior_size
        end
    end
    
    return alpha
end

"""
    extract_pin_config(state::UInt32, pins::Vector{Int})

Extract the configuration of pin vertices from a full state.

# Arguments
- `state`: Full vertex configuration as bit mask
- `pins`: Pin vertex indices (1-based)

# Returns
- `UInt32`: Pin configuration as bit mask (pin i → bit i-1)
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
    verify_gadget_via_alpha_tensor(
        graph::SimpleGraph{Int},
        pins::Vector{Int},
        target_states::Vector{UInt32}
    ) -> Union{Nothing, Tuple{Vector{Float64}, Int}}

Verify if a graph can serve as a gadget using α-tensor equivalence.

Returns uniform weights (all 1.0) if the gadget is valid, or nothing if invalid.
Also returns the MIS overhead constant c.

# Returns
- `nothing` if gadget is invalid
- `(weights::Vector{Float64}, overhead::Int)` if valid
"""
function verify_gadget_via_alpha_tensor(
    graph::SimpleGraph{Int},
    pins::Vector{Int},
    target_states::Vector{UInt32}
)
    # 1. Infer pattern's α-tensor from target states
    α_pattern = infer_pattern_alpha(target_states, pins)
    
    # 2. Compute gadget's α-tensor
    α_gadget = compute_reduced_alpha_tensor(graph, pins)
    
    # 3. Check equivalence
    is_equiv, overhead = check_alpha_equivalence(α_pattern, α_gadget)
    
    if !is_equiv
        return nothing
    end
    
    # Valid gadget - return uniform weights
    n_vertices = Graphs.nv(graph)
    weights = ones(Float64, n_vertices)
    
    return (weights, overhead)
end

