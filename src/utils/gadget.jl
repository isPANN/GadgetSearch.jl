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
    analyze_gadget(gadget::Gadget; model::Type{<:EnergyModel}=RydbergModel)

Compute the maximum-energy states of a gadget and return structured data for
programmatic consumers such as the visual editor.
"""
function analyze_gadget(gadget::Gadget; model::Type{<:EnergyModel}=RydbergModel)
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
    state_count > 0 || error("No states found.")

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

    model_name = model === RydbergModel ? "Rydberg (MIS)" : "QUBO (Full)"
    ground_states = map(max_indices) do idx
        config = states[idx]
        (
            state_index=idx,
            configuration=Int[((config >> (v - 1)) & 0x1) for v in 1:num_vertices],
            pins=Int[((config >> (p - 1)) & 0x1) for p in pins],
        )
    end

    return (
        model=model_name,
        state_count=state_count,
        max_energy=max_energy,
        ground_states=ground_states,
    )
end

"""
    check_gadget(gadget::Gadget; _return_info::Bool=false, model::Type{<:EnergyModel}=RydbergModel)

Validate a `Gadget` by computing energies for its state space and
reporting the ground state configurations on pins.
"""
function check_gadget(gadget::Gadget; _return_info::Bool=false, model::Type{<:EnergyModel}=RydbergModel)
    report = analyze_gadget(gadget; model=model)
    lines = [
        "Model: $(report.model)",
        "Max energy value: $(report.max_energy)",
        "Ground states (max energy):",
    ]
    for state in report.ground_states
        push!(lines, "  State index=$(state.state_index), pins=$(state.pins)")
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
