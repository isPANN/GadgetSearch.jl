# Gadget utility functions
# The Gadget struct is defined in src/core/search.jl

function save_results_to_json(results::Vector{<:Gadget}, file_path::String)
    # Convert each result to a serializable dictionary
    json_results = map(results) do res
        base_dict = Dict(
            "ground_states" => [vec(res.ground_states[i, :]) for i in 1:size(res.ground_states, 1)],
            "pins" => res.pins,
            "graph" => Dict(
                "nodes" => [Dict("id" => i, "weight" => res.vertex_weights[i]) for i in 1:length(res.vertex_weights)],
                "edges" => [Dict("source" => src(e), "target" => dst(e)) for e in Graphs.edges(res.graph)],
            )
        )

        # Add edge weights for QUBO gadgets
        if !isempty(res.edge_weights)
            base_dict["edge_weights"] = [
                Dict("source" => res.edge_list[i][1], "target" => res.edge_list[i][2], "weight" => res.edge_weights[i])
                for i in 1:length(res.edge_weights)
            ]
        end

        # Add position data if available
        if !isnothing(res.pos)
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

"""
    check_gadget(gadget::Gadget; _return_info::Bool=false, model::Type{<:EnergyModel}=RydbergModel)

Validate a `Gadget` by computing energies for its state space and
reporting the ground state configurations on pins.

# Arguments
- `gadget::Gadget`: The gadget to check.

# Keyword Arguments
- `_return_info::Bool=false`: If `true`, return a string; otherwise log with `@info`.
- `model::Type{<:EnergyModel}=RydbergModel`: Energy model to use for validation.

# Returns
- When `_return_info` is `true`, returns a formatted `String`; otherwise returns `nothing`.
"""
function check_gadget(gadget::Gadget; _return_info::Bool=false, model::Type{<:EnergyModel}=RydbergModel)
    g = gadget.graph
    vertex_weights = gadget.vertex_weights
    edge_weights = gadget.edge_weights
    edge_list = gadget.edge_list
    pins = gadget.pins

    num_vertices = nv(g)
    length(vertex_weights) == num_vertices || throw(ArgumentError("length(vertex_weights)=$(length(vertex_weights)) must equal nv(g)=$num_vertices"))
    all(1 .<= pins .<= num_vertices) || throw(ArgumentError("pins must be within 1:nv(g)"))

    # Get state space based on model
    states, state_count = get_state_space(model, g)
    state_count > 0 || return _return_info ? "No states found." : (@info "No states found."; nothing)

    # Helper to compute energy of a state
    function _energy_of_config(config::Unsigned)
        # Vertex energy
        vertex_energy = zero(float(eltype(vertex_weights)))
        for v in 1:num_vertices
            if ((config >> (v - 1)) & 0x1) == 1
                vertex_energy += vertex_weights[v]
            end
        end
        
        # Edge energy (for QUBO)
        edge_energy = zero(float(eltype(vertex_weights)))
        for (i, (u, v)) in enumerate(edge_list)
            if ((config >> (u - 1)) & 0x1) == 1 && ((config >> (v - 1)) & 0x1) == 1
                edge_energy += edge_weights[i]
            end
        end
        
        return vertex_energy + edge_energy
    end

    # Compute energies
    WeightFloat = float(eltype(vertex_weights))
    energy_values = Vector{WeightFloat}(undef, state_count)
    @inbounds for i in 1:state_count
        energy_values[i] = _energy_of_config(states[i])
    end

    # Find maxima (for MIS we maximize, for QUBO we might want to check ground states differently)
    max_energy = maximum(energy_values)
    # Use approximate equality for floating point comparison (tolerance 1e-6)
    max_indices = findall(e -> abs(e - max_energy) < 1e-6, energy_values)

    # Format report
    model_name = model === RydbergModel ? "Rydberg (MIS)" : model === RydbergUnweightedModel ? "RydbergUnweighted (α-tensor MIS)" : "QUBO (Full)"
    lines = ["Model: $model_name", "Max energy value: $(max_energy)", "Ground states (max energy):"]
    for idx in max_indices
        config = states[idx]
        pin_values = @inbounds Int[((config >> (p - 1)) & 0x1) for p in pins]
        push!(lines, "  State index=$(idx), pins=$(pin_values)")
    end

    msg = join(lines, "\n")
    return _return_info ? msg : (@info msg; nothing)
end

"""
    check_gadget_rydberg(gadget::Gadget; _return_info::Bool=false)

Check gadget using Rydberg (MIS) model.
"""
check_gadget_rydberg(gadget::Gadget; kwargs...) = check_gadget(gadget; model=RydbergModel, kwargs...)

"""
    check_gadget_qubo(gadget::Gadget; _return_info::Bool=false)

Check gadget using QUBO (full state space) model.
"""
check_gadget_qubo(gadget::Gadget; kwargs...) = check_gadget(gadget; model=QUBOModel, kwargs...)

"""
    check_gadget_unweighted(gadget::Gadget; _return_info::Bool=false)

Check gadget using unweighted Rydberg (α-tensor MIS) model with uniform weights.
Ground states are the Maximum Independent Sets (largest cardinality MIS).
Also computes and displays the reduced α-tensor for the gadget.
"""
function check_gadget_unweighted(gadget::Gadget; _return_info::Bool=false)
    # Standard check with unweighted model  
    result = check_gadget(gadget; model=RydbergUnweightedModel, _return_info=true)
    
    # Additionally compute α-tensor info
    g = gadget.graph
    pins = gadget.pins
    α = compute_reduced_alpha_tensor(g, pins)
    n_pins = length(pins)
    ground_configs, max_total = find_ground_configs(α, n_pins)
    
    alpha_lines = ["\n--- α-Tensor Analysis ---", "Max total IS size: $max_total"]
    push!(alpha_lines, "Ground configs (boundary configs with max total):")
    for config in sort(collect(ground_configs))
        pin_vals = Int[((config >> (i-1)) & 0x1) for i in 1:n_pins]
        push!(alpha_lines, "  config=$(pin_vals), α=$(α[config]), total=$(α[config] + count_ones(config))")
    end
    push!(alpha_lines, "All α-tensor entries:")
    for config in UInt32(0):UInt32(2^n_pins - 1)
        α_val = α[config]
        pin_vals = Int[((config >> (i-1)) & 0x1) for i in 1:n_pins]
        status = α_val == INFEASIBLE_ALPHA ? "INFEASIBLE" : "α=$(α_val), total=$(α_val + count_ones(config))"
        push!(alpha_lines, "  config=$(pin_vals): $status")
    end
    
    msg = result * "\n" * join(alpha_lines, "\n")
    return _return_info ? msg : (@info msg; nothing)
end
